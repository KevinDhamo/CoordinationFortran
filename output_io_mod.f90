!> @brief Module for writing coordination analysis results to output files
module output_io_mod
    use types_mod
    use config_mod
    use error_mod
    use coordination_mod
    use cell_list_mod, only: num_cell_resets, num_cell_updates, num_selective_updates, &
                            get_cell_statistics, total_empty_cells, &
                            max_atoms_in_cell, min_atoms_in_cell, &
                            avg_atoms_per_cell, n_cells_x, n_cells_y, n_cells_z
    use atom_selection_mod, only: is_atom_selected, selection_active
    use omp_lib
    implicit none
    private

    !===============================================================================
    ! Public procedures
    !===============================================================================
    public :: init_output
    public :: write_output_frame
    public :: write_final_report
    public :: cleanup_output
    
    !===============================================================================
    ! Output control parameters
    !===============================================================================
    !> Option to disable output file writing (e.g., for benchmarking)
    logical, public :: disable_output = .false.

    !===============================================================================
    ! Module variables
    !===============================================================================
    !> File unit for output file
    integer :: output_unit = -1
    
    !> Flag indicating if output file is open
    logical :: output_file_open = .false.

contains
    !> @brief Initialize output file
    !>
    !> Opens the output file and writes the initial header with column descriptions.
    !>
    !> @param[in] filename Path to output file
    subroutine init_output(filename)
        character(len=*), intent(in) :: filename
        integer :: io_stat, i
        
        ! If output is disabled, do nothing
        if (disable_output) return
        
        ! Open output file
        open(newunit=output_unit, file=filename, status='replace', &
             action='write', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open output file: " // trim(filename), ERR_FILE_IO)
        end if
        output_file_open = .true.
        
        ! Write file header with column explanations
        write(output_unit,'(A)') "# Coordination Number Analysis"
        write(output_unit,'(A)') "# Column format:"
        write(output_unit,'(A)') "#   1. Atom ID"
        write(output_unit,'(A)') "#   2. Atom Type"
        
        ! Write pair column descriptions
        do i = 1, n_pairs
            write(output_unit,'(A,I0,A,I0,A,I0)') "#   ", i+2, &
                ". CN(Type ", pairs(i)%type1_id, "-", pairs(i)%type2_id, ")"
        end do
        
        ! Add note if atom selection is active
        if (selection_active) then
            write(output_unit,'(A)') "# Note: Only selected atoms are included in this output"
        end if
        
        ! Add note if neighbor info is shown
        if (show_neighbors) then
            write(output_unit,'(A)') "# Note: Neighbor atom IDs are shown in parentheses after coordination numbers"
        end if
        
        ! Add note about calculation method
        if (use_cell_list) then
            write(output_unit,'(A)') "# Note: Using Cell-List method for coordination calculation"
        else
            write(output_unit,'(A)') "# Note: Using direct O(n²) method for coordination calculation"
        end if
        
        write(output_unit,'(A)') "# Frame data follows"
    end subroutine init_output

    !> @brief Write coordination data for a trajectory frame
    !>
    !> Writes the coordination numbers for all atoms in the current frame
    !> to the output file. If neighbor tracking is enabled, also writes
    !> the IDs of neighboring atoms in a consistent, parsable format.
    !>
    !> @param[in] frame_number Current frame number
    !> @param[in] elements Array of element names for each atom
    !> @param[in] atom_types Array of atom type indices
    !> @param[in] coord_numbers_local Coordination numbers matrix (atom_index, pair_index)
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] atom_info Array of atom type information
    subroutine write_output_frame(frame_number, elements, atom_types, coord_numbers_local, n_atoms, atom_info)
        integer, intent(in) :: frame_number
        character(len=*), intent(in) :: elements(:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: coord_numbers_local(:,:)
        integer, intent(in) :: n_atoms
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: i, j, k, neighbor_id
        character(len=32) :: fmt_str
        character(len=2048) :: buffer
        character(len=20) :: temp_str
        
        ! Skip if output is disabled
        if (disable_output .or. .not. output_file_open) return
        
        ! Write frame header
        write(output_unit,'(A,I0)') "# Frame: ", frame_number
        
        ! Handle output format based on whether we show neighbors or not
        if (.not. show_neighbors) then
            ! Standard output format without neighbors
            write(fmt_str,'(A,I0,A)') '(I6,A15,', n_pairs, 'I10)'
            
            ! Write atom data with type names - only for selected atoms
            do i = 1, n_atoms
                ! Skip atoms that are not selected
                if (.not. is_atom_selected(i)) cycle
                
                if (atom_types(i) >= 1 .and. atom_types(i) <= size(atom_info)) then
                    write(output_unit, fmt_str) i, trim(atom_info(atom_types(i))%name), &
                                                (coord_numbers_local(i,j), j=1,n_pairs)
                else
                    write(output_unit, fmt_str) i, "Unknown", (coord_numbers_local(i,j), j=1,n_pairs)
                end if
            end do
        else
            ! Enhanced output format with neighbor atom IDs - only for selected atoms
            do i = 1, n_atoms
                ! Skip atoms that are not selected
                if (.not. is_atom_selected(i)) cycle
                
                ! Write atom ID and type
                if (atom_types(i) >= 1 .and. atom_types(i) <= size(atom_info)) then
                    write(buffer, '(I6,A15)') i, trim(atom_info(atom_types(i))%name)
                else
                    write(buffer, '(I6,A15)') i, "Unknown"
                end if
                
                ! Add coordination numbers with neighbor lists
                do j = 1, n_pairs
                    ! Add coordination number with consistent spacing
                    write(temp_str, '(I6)') coord_numbers_local(i,j)
                    buffer = trim(buffer) // " " // adjustl(temp_str)
                    
                    ! Add neighbor atom IDs if available
                    if (coord_numbers_local(i,j) > 0) then
                        ! Always use the same format with no spaces to make parsing easier
                        buffer = trim(buffer) // "("
                        
                        ! Add up to MAX_NEIGHBORS neighbor IDs
                        do k = 1, min(coord_numbers_local(i,j), size(coord_neighbors, 3))
                            neighbor_id = coord_neighbors(i, j, k)
                            if (neighbor_id > 0) then
                                if (k > 1) buffer = trim(buffer) // ","
                                write(temp_str, '(I0)') neighbor_id
                                buffer = trim(buffer) // trim(adjustl(temp_str))
                            else
                                ! We've reached the end of valid neighbors
                                exit
                            end if
                        end do
                        buffer = trim(buffer) // ")"
                    end if
                end do
                
                ! Write the complete line to output
                write(output_unit, '(A)') trim(buffer)
            end do
        end if
    end subroutine write_output_frame

    !> @brief Write final analysis report with statistics
    !>
    !> Generates a comprehensive report with coordination statistics,
    !> cell list performance metrics, and system composition information.
    !>
    !> @param[in] n_frames Number of frames processed
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] total_time Total processing time in seconds
    !> @param[in] atom_info Array of atom type information
    !> @param[in] n_types Number of atom types
    !> @param[in] coord_numbers_local Final coordination numbers matrix
    !> @param[in] atom_types_local Array of atom type indices
    subroutine write_final_report(n_frames, n_atoms, total_time, atom_info, &
                                 n_types, coord_numbers_local, atom_types_local)
        integer, intent(in) :: n_frames, n_atoms, n_types
        real(dp), intent(in) :: total_time
        type(atom_type_info), intent(in) :: atom_info(:)
        integer, intent(in) :: coord_numbers_local(:,:)
        integer, intent(in) :: atom_types_local(:)
        
        real(dp), allocatable :: type_avg_coord(:,:)
        integer, allocatable :: type_counts(:)
        integer :: i, j, k, alloc_stat
        real(dp) :: cell_stats(4)
        real(dp) :: empty_cell_percentage
        logical, allocatable :: pair_in_setup(:)
        real(dp), allocatable :: avg_pair_values(:)
        integer :: type1, type2   ! Add these declarations
        
        ! Allocate arrays for statistics with error checking
        allocate(type_avg_coord(n_types, n_pairs), type_counts(n_types), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate statistics arrays", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize arrays
        type_avg_coord = 0.0_dp
        type_counts = 0
        
        ! Calculate statistics safely - only for selected atoms
        !$OMP PARALLEL DO REDUCTION(+:type_counts,type_avg_coord) PRIVATE(i,j) SCHEDULE(guided)
        do i = 1, n_atoms
            ! Skip atoms that were not selected
            if (.not. is_atom_selected(i)) cycle
            
            if (atom_types_local(i) > 0 .and. atom_types_local(i) <= n_types) then
                if (atom_info(atom_types_local(i))%include) then
                    type_counts(atom_types_local(i)) = type_counts(atom_types_local(i)) + 1
                    do j = 1, n_pairs
                        type_avg_coord(atom_types_local(i),j) = &
                            type_avg_coord(atom_types_local(i),j) + coord_numbers_local(i,j)
                    end do
                end if
            end if
        end do
        !$OMP END PARALLEL DO
        
        ! Calculate averages
        do i = 1, n_types
            if (type_counts(i) > 0) then
                type_avg_coord(i,:) = type_avg_coord(i,:) / type_counts(i)
            end if
        end do
        
        ! Only report cell list statistics if using cell list method
        if (use_cell_list) then
            ! Get cell list statistics safely
            cell_stats = get_cell_statistics()
            
            ! Calculate percentage of empty cells
            empty_cell_percentage = (total_empty_cells * 100.0) / &
                                    max(1, n_cells_x * n_cells_y * n_cells_z)
            
            !===============================================================================
            ! Cell list statistics
            !===============================================================================
            write(*,'(A)') " CELL LIST STATISTICS"
            write(*,'(A)') " -------------------"
            write(*,'(A,I0)') " Complete cell list rebuilds:      ", num_cell_resets
            write(*,'(A,I0)') " Cell list maintenance updates:    ", num_cell_updates
            write(*,'(A,I0)') " Selective cell rebuilds:          ", num_selective_updates
            write(*,'(A,I0,A,F5.1,A)') " Average empty cells:              ", total_empty_cells, &
                                     " (", empty_cell_percentage, "%)"
            write(*,'(A,F10.2)') " Average atoms per non-empty cell: ", avg_atoms_per_cell
            write(*,'(A,I0,A,I0)') " Min/Max atoms per cell:           ", min_atoms_in_cell, "/", max_atoms_in_cell
            write(*,*)
        end if
        
        ! Allocate arrays for tracking which pairs were defined in setup
        allocate(pair_in_setup(n_pairs), avg_pair_values(n_pairs), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate pair tracking arrays", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize pair arrays
        pair_in_setup = .false.
        avg_pair_values = 0.0_dp
        
        ! Mark pairs that were defined in setup (not using default cutoff)
        do j = 1, n_pairs
            type1 = pairs(j)%type1_id
            type2 = pairs(j)%type2_id
            
            ! Skip pairs with empty type names or no atoms
            if (len_trim(atom_info(type1)%name) == 0 .or. len_trim(atom_info(type2)%name) == 0) cycle
            if (type_counts(type1) == 0 .or. type_counts(type2) == 0) cycle
            
            ! Consider pair as defined in setup if it doesn't use default cutoff or is self-interaction
            if (pairs(j)%cutoff /= DEFAULT_CUTOFF .or. type1 == type2) then
                pair_in_setup(j) = .true.
            end if
        end do
        
        !===============================================================================
        ! Coordination statistics - simplified to show only unique pair types
        !===============================================================================
        write(*,'(A)') " COORDINATION STATISTICS"
        write(*,'(A)') " ----------------------"
        
        ! Display atom counts by type first
        write(*,'(A)') " Atom counts by type:"
        k = 0  ! Counter for valid types
        do i = 1, n_types
            if (atom_info(i)%include .and. type_counts(i) > 0 .and. len_trim(atom_info(i)%name) > 0) then
                write(*,'(A,I0,A,A,A,I0,A)') "   Type ", i, " (", trim(atom_info(i)%name), "): ", &
                    type_counts(i), " atoms"
                k = k + 1
            end if
        end do
        
        if (k == 0) then
            write(*,'(A)') "   No valid atom types found"
        end if
        write(*,*)
        
        ! Calculate average coordination values for each pair type
        do j = 1, n_pairs
            type1 = pairs(j)%type1_id
            type2 = pairs(j)%type2_id
            
            ! Only process valid pairs
            if (pair_in_setup(j) .and. type_counts(type1) > 0 .and. type_counts(type2) > 0) then
                if (type1 == type2) then
                    ! For same-type pairs (like Na-Na), just use the direct average
                    avg_pair_values(j) = type_avg_coord(type1, j)
                else
                    ! For different types (like Na-Cl), calculate weighted average
                    avg_pair_values(j) = (type_avg_coord(type1, j) * type_counts(type1) + &
                                         type_avg_coord(type2, j) * type_counts(type2)) / &
                                         (type_counts(type1) + type_counts(type2))
                end if
            end if
        end do
        
        ! Display average coordination numbers
        write(*,'(A)') " Average coordination numbers by pair type:"
        k = 0  ! Counter for valid pairs
        do j = 1, n_pairs
            type1 = pairs(j)%type1_id
            type2 = pairs(j)%type2_id
            
            ! Only show pairs that were defined in setup and have valid types with atoms
            if (pair_in_setup(j) .and. &
                len_trim(atom_info(type1)%name) > 0 .and. &
                len_trim(atom_info(type2)%name) > 0 .and. &
                type_counts(type1) > 0 .and. type_counts(type2) > 0) then
                
                write(*,'(A,A,A,I0,A,A,A,I0,A,F10.3)') &
                    "   ", trim(atom_info(type1)%name), "-", &
                    type1, " to ", &
                    trim(atom_info(type2)%name), "-", &
                    type2, ":     ", avg_pair_values(j)
                k = k + 1
            end if
        end do
        
        if (k == 0) then
            write(*,'(A)') "   No valid pair interactions defined"
        end if
        write(*,*)
        
        ! Clean up arrays
        if (allocated(type_avg_coord)) deallocate(type_avg_coord)
        if (allocated(type_counts)) deallocate(type_counts)
        if (allocated(pair_in_setup)) deallocate(pair_in_setup)
        if (allocated(avg_pair_values)) deallocate(avg_pair_values)
    end subroutine write_final_report

    !> @brief Clean up output resources
    !>
    !> Closes the output file if it's open.
    subroutine cleanup_output()
        if (output_file_open) then
            close(output_unit)
            output_file_open = .false.
        end if
    end subroutine cleanup_output

end module output_io_mod
