!> @brief Module for writing coordination analysis results to output files
module output_io_mod
    use types_mod
    use config_mod
    use error_mod
    use coordination_mod, only: pairs, n_pairs, coord_neighbors, show_neighbors, verlet_stats
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
        if (use_verlet_lists) then
            write(output_unit,'(A)') "# Note: Using hybrid Cell-List/Verlet-List method for coordination calculation"
        else
            write(output_unit,'(A)') "# Note: Using traditional Cell-List method for coordination calculation"
        end if
        
        write(output_unit,'(A)') "# Frame data follows"
    end subroutine init_output

    !> @brief Write coordination data for a trajectory frame
    !>
    !> Writes the coordination numbers for all atoms in the current frame
    !> to the output file. If neighbor tracking is enabled, also writes
    !> the IDs of neighboring atoms.
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
                    ! Add coordination number - first one gets a space
                    if (j == 1) then
                        write(temp_str, '(I5)') coord_numbers_local(i,j)
                        buffer = trim(buffer) // " " // trim(adjustl(temp_str))
                    else
                        write(temp_str, '(I5)') coord_numbers_local(i,j)
                        buffer = trim(buffer) // "   " // trim(adjustl(temp_str))
                    end if
                    
                    ! Add neighbor atom IDs if available
                    if (coord_numbers_local(i,j) > 0) then
                        ! First pair gets a space, others no space before parenthesis
                        if (j == 1) then
                            buffer = trim(buffer) // " ("
                        else
                            buffer = trim(buffer) // "("
                        end if
                        
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
        integer :: i, j, alloc_stat
        real(dp) :: cell_stats(4)
        real(dp) :: empty_cell_percentage
        
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
        
        !===============================================================================
        ! Coordination statistics
        !===============================================================================
        write(*,'(A)') " COORDINATION STATISTICS"
        write(*,'(A)') " ----------------------"
        do i = 1, n_types
            if (atom_info(i)%include) then
                write(*,'(A,I0,A,A,A,I0,A)') " Type ", i, " (", trim(atom_info(i)%name), ") - ", &
                    type_counts(i), " atoms:"
                    
                do j = 1, n_pairs
                    if (pairs(j)%type1_id == i .or. pairs(j)%type2_id == i) then
                        ! Format according to ideal output format
                        write(*,'(A,A,A,I0,A,A,A,I0,A,F10.3)') &
                            "   With ", trim(atom_info(pairs(j)%type1_id)%name), "-", &
                            pairs(j)%type1_id, " to ", &
                            trim(atom_info(pairs(j)%type2_id)%name), "-", &
                            pairs(j)%type2_id, ":     ", type_avg_coord(i,j)
                    end if
                end do
                write(*,*)
            end if
        end do
        
        ! Clean up arrays
        if (allocated(type_avg_coord)) deallocate(type_avg_coord)
        if (allocated(type_counts)) deallocate(type_counts)
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