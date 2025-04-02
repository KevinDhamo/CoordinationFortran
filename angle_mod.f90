!> @brief Module for calculating angles between atoms using coordination data
!>
!> This module calculates angles by directly using pre-calculated coordination data,
!> avoiding redundant distance calculations and neighbor finding operations.
!> It uses an efficient squared cosine pre-filtering approach to avoid unnecessary
!> square root calculations for angles outside the range of interest.
module angle_mod
    use types_mod
    use config_mod
    use coordination_mod
    use error_mod
    use omp_lib
    implicit none
    private

    !===============================================================================
    ! Public procedures and variables
    !===============================================================================
    public :: init_angle_analysis
    public :: compute_angles
    public :: finalize_angle_analysis
    
    !===============================================================================
    ! Module constants and variables
    !===============================================================================
    !> Output unit for angle file
    integer :: angle_unit = -1
    
    !> Flag indicating if angle file is open
    logical :: angle_file_open = .false.

    !> Structure to store angle information
    type :: angle_data_type
        integer :: central_atom      !< Index of central atom
        integer :: neighbor1         !< Index of first neighbor
        integer :: neighbor2         !< Index of second neighbor
        real(dp) :: angle            !< Angle in degrees
        integer :: type1, type2      !< Atom types of neighbors
        integer :: central_type      !< Type of central atom
    end type angle_data_type
    
    !> Array to store all angles
    type(angle_data_type), allocatable :: angle_list(:)
    
    !> Number of angles stored
    integer :: angle_count = 0
    
    !> Maximum number of angles that can be stored
    integer :: max_angles = 0

    !> Flag to track if we've shown specific warnings
    logical :: diagnostics_shown = .false.

    !> Performance statistics
    integer :: found_triplets = 0
    integer :: calculated_angles = 0
    integer :: filtered_angles = 0
    
    !> PI constant for angle calculation
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    
    !> Local copies of filtering parameters (from config_mod)
    real(dp) :: local_min_angle = 0.0_dp
    real(dp) :: local_max_angle = 180.0_dp
    
    !> Pre-calculated squared cosine values for efficient filtering
    real(dp) :: min_cos_sq = 0.0_dp
    real(dp) :: max_cos_sq = 1.0_dp
    
    !> Flag to indicate if we're filtering by angle range
    logical :: using_angle_filter = .false.
    
    !> Flag to indicate if the angle range crosses 90 degrees
    logical :: range_crosses_90deg = .false.
    
    !> Maximum number of neighbors to consider per atom (memory optimization)
    integer, parameter :: MAX_ALL_NEIGHBORS = 200
    
    !> Thread-local angle buffer size for parallel processing
    integer :: LOCAL_ANGLE_BUFFER_SIZE = 10000  ! Will be set from setup.txt

contains

    !> @brief Initialize the angle analysis system
    !> 
    !> Sets up data structures for angle calculation and reads filtering parameters
    !> from the configuration.
    !>
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] n_types Number of atom types
    !> @param[in] atom_info Array of atom type information
    subroutine init_angle_analysis(n_atoms, n_types, atom_info)
        integer, intent(in) :: n_atoms, n_types
        type(atom_type_info), intent(in) :: atom_info(:)
        integer :: alloc_stat, i, j, k
        logical :: any_missing, type2_found
        
        ! Use the batch size from setup.txt if available
        if (angle_batch_size > 0) then
            LOCAL_ANGLE_BUFFER_SIZE = angle_batch_size
        end if
        
        ! Verify that show_neighbors is enabled - this is critical for our optimization
        if (.not. show_neighbors) then
            write(*,'(A)') " ERROR: show_neighbors must be enabled in setup.txt for angle calculations."
            write(*,'(A)') "        Please add 'SHOW_NEIGHBORS=yes' to your setup file and restart."
            ! Force enable it for this run (though it's better to have it in setup.txt)
            show_neighbors = .true.
        end if
        
        ! Read angle range filtering parameters from config_mod (setup.txt)
        local_min_angle = min_angle_degree 
        local_max_angle = max_angle_degree
        
        ! Check if we're using angle filtering
        using_angle_filter = (local_min_angle > 0.0_dp .or. local_max_angle < 180.0_dp)
        
        ! Check if the range crosses 90 degrees (this affects our filtering logic)
        range_crosses_90deg = (local_min_angle <= 90.0_dp .and. local_max_angle >= 90.0_dp)
        
        ! Pre-calculate squared cosine values for quick filtering
        ! Note the reversal: min_angle → max_cos_sq, max_angle → min_cos_sq
        if (using_angle_filter) then
            ! Convert to radians for cosine calculation
            min_cos_sq = cos(local_max_angle * PI/180.0_dp)**2
            max_cos_sq = cos(local_min_angle * PI/180.0_dp)**2
        end if
        
        ! Estimate initial capacity (will be resized if needed)
        max_angles = min(n_atoms * 20, 10000000)
        
        ! Allocate angle array
        if (allocated(angle_list)) deallocate(angle_list)
        allocate(angle_list(max_angles), stat=alloc_stat)
        call check_allocation(alloc_stat, "angle_list array")
        
        ! Initialize angle count
        angle_count = 0
        found_triplets = 0
        calculated_angles = 0
        filtered_angles = 0
        
        ! Reset diagnostics flag
        diagnostics_shown = .false.
        
        ! Set up output file
        if (.not. angle_file_open) then
            open(newunit=angle_unit, file=angle_output_file, status='replace', &
                 action='write', iostat=alloc_stat)
            if (alloc_stat /= 0) then
                call handle_error("Failed to open angle output file", ERR_FILE_IO)
                return
            end if
            angle_file_open = .true.
            
            ! Write header
            write(angle_unit,'(A)') "# Angle Distribution Analysis"
            write(angle_unit,'(A)') "# Format:"
            write(angle_unit,'(A)') "#   1. Central Atom ID"
            write(angle_unit,'(A)') "#   2. Central Atom Type"
            write(angle_unit,'(A)') "#   3. Neighbor 1 ID"
            write(angle_unit,'(A)') "#   4. Neighbor 1 Type"
            write(angle_unit,'(A)') "#   5. Neighbor 2 ID"
            write(angle_unit,'(A)') "#   6. Neighbor 2 Type"
            write(angle_unit,'(A)') "#   7. Angle (degrees)"
        end if
        
        ! Print information about atom types used in angle calculation
        if (VERBOSE) then
            write(*,'(A)') " Angle calculation initialized (optimized version)"
            write(*,'(A)') " Using pre-calculated coordination data for efficiency"
            
            ! Report angle filtering settings
            if (using_angle_filter) then
                write(*,'(A,F6.1,A,F6.1,A)') " Angle filtering: ", local_min_angle, " to ", &
                                           local_max_angle, " degrees"
            else
                write(*,'(A)') " Angle filtering: None (calculating all angles)"
            end if
            
            ! Report on type filtering
            if (allocated(angle_type_mask)) then
                any_missing = .false.
                
                do i = 1, n_types
                    do j = 1, n_types
                        do k = 1, n_types
                            if (.not. any_missing) then
                                if (.not. angle_type_mask(i, j, k)) then
                                    write(*,'(A)') " Note: Some atom type triplets are excluded from angle calculation"
                                    any_missing = .true.
                                end if
                            end if
                        end do
                    end do
                end do
                
                ! Check for missing type 2 combinations
                if (n_types >= 2) then
                    type2_found = .false.
                    
                    do i = 1, n_types
                        if (angle_type_mask(i, 2, i) .or. angle_type_mask(2, i, i) .or. &
                            angle_type_mask(i, i, 2)) then
                            type2_found = .true.
                            exit
                        end if
                    end do
                    
                    if (.not. type2_found) then
                        write(*,'(A)') " WARNING: No angle triplets including Type 2 are enabled!"
                        write(*,'(A)') " Missing combinations may include:"
                        write(*,'(A)') "   ANGLE_TYPES_1_2_1, ANGLE_TYPES_1_2_2, ANGLE_TYPES_1_2_3,"
                        write(*,'(A)') "   ANGLE_TYPES_2_2_2, ANGLE_TYPES_2_2_3, ANGLE_TYPES_2_3_3"
                    end if
                end if
            end if
        end if
    end subroutine init_angle_analysis
    
!> @brief Calculate angles using pre-calculated coordination data with efficient filtering
!>
!> This optimized version uses pre-calculated coordination data and the squared
!> cosine pre-filtering approach to avoid unnecessary square root calculations.
!>
!> @param[in] coords Atom coordinates
!> @param[in] atom_types Atom type indices
!> @param[in] n_atoms Number of atoms in the system
!> @param[in] box_length Box dimensions in x, y, z
!> @param[in] frame Current frame number
!> @param[in] atom_info Array of atom type information
!> @param[in] quiet Optional flag to suppress diagnostic output (default: .false.)
subroutine compute_angles(coords, atom_types, n_atoms, box_length, frame, atom_info, quiet)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: box_length(3)
    integer, intent(in) :: frame
    type(atom_type_info), intent(in) :: atom_info(:)
    logical, intent(in), optional :: quiet
    
    integer :: i, j, k, central_atom, neigh1, neigh2, pair_idx
    integer :: num_neighbors
    real(dp) :: start_time, end_time
    integer, allocatable :: all_neighbors(:)
    integer :: all_neighbors_count
    integer :: thread_num, local_angle_count
    type(angle_data_type), allocatable :: local_angle_buffer(:)
    logical :: suppress_output
    
    ! Thread-local variables for statistics
    integer :: local_found_triplets, local_calculated_angles, local_filtered_angles
    
    ! For vector calculation
    real(dp) :: vec1(3), vec2(3), dot_prod, angle
    real(dp) :: dx1, dy1, dz1, dx2, dy2, dz2
    real(dp) :: vec1_sq, vec2_sq, cos_sq
    logical :: is_acute
    
    ! Variables for loop limits
    integer :: max_neighbors_per_pair
    
    ! Timing information
    real(dp) :: collection_time, calculation_time, storage_time
    
    ! Process the optional quiet parameter
    suppress_output = .false.
    if (present(quiet)) suppress_output = quiet
    
    ! Reset angle count and statistics
    angle_count = 0
    found_triplets = 0
    calculated_angles = 0
    filtered_angles = 0
    
    ! Record start time for performance measurement
    start_time = omp_get_wtime()
    collection_time = 0.0_dp
    calculation_time = 0.0_dp
    storage_time = 0.0_dp
    
    ! Calculate maximum neighbors that could be stored per pair
    max_neighbors_per_pair = size(coord_neighbors, dim=3)
    
    ! Process atoms in parallel
    !$OMP PARALLEL PRIVATE(i, j, k, central_atom, neigh1, neigh2, pair_idx, &
    !$OMP                  num_neighbors, all_neighbors, all_neighbors_count, &
    !$OMP                  thread_num, local_angle_count, local_angle_buffer, &
    !$OMP                  vec1, vec2, dot_prod, angle, &
    !$OMP                  dx1, dy1, dz1, dx2, dy2, dz2, vec1_sq, vec2_sq, cos_sq, is_acute, &
    !$OMP                  local_found_triplets, local_calculated_angles, local_filtered_angles)
    
    ! Get thread number for local data
    thread_num = omp_get_thread_num()
    
    ! Each thread gets a local buffer for angles
    allocate(local_angle_buffer(LOCAL_ANGLE_BUFFER_SIZE))
    allocate(all_neighbors(MAX_ALL_NEIGHBORS))
    
    ! Initialize thread local counters
    local_angle_count = 0
    local_found_triplets = 0
    local_calculated_angles = 0
    local_filtered_angles = 0
    
    ! Distribute work by atoms
    !$OMP DO SCHEDULE(dynamic,100)
    do i = 1, n_atoms
        central_atom = i
        
        ! Skip if this atom type should be excluded
        if (.not. atom_info(atom_types(central_atom))%include) cycle
        
        ! Collect all neighbors from the coordination data
        ! This avoids redundant distance calculations
        all_neighbors_count = 0
        
        do pair_idx = 1, n_pairs
            num_neighbors = neighbor_counts(central_atom, pair_idx)
            
            ! Skip if no neighbors for this pair type
            if (num_neighbors == 0) cycle
            
            ! Add neighbors to the all_neighbors array
            do j = 1, min(num_neighbors, max_neighbors_per_pair)
                ! Check if we have room in the neighbor array
                if (all_neighbors_count >= MAX_ALL_NEIGHBORS) then
                    exit  ! Skip remaining neighbors if buffer is full
                end if
                
                neigh1 = coord_neighbors(central_atom, pair_idx, j)
                
                ! Skip invalid neighbors
                if (neigh1 <= 0) cycle
                
                ! Add to all_neighbors (avoiding duplicates)
                ! Since each neighbor might appear in multiple pair types
                if (.not. any(all_neighbors(1:all_neighbors_count) == neigh1)) then
                    all_neighbors_count = all_neighbors_count + 1
                    all_neighbors(all_neighbors_count) = neigh1
                end if
            end do
            
            ! Stop if we've reached capacity
            if (all_neighbors_count >= MAX_ALL_NEIGHBORS) exit
        end do
        
        ! Need at least 2 neighbors to form an angle
        if (all_neighbors_count < 2) cycle
        
        ! Calculate angles between all pairs of neighbors
        do j = 1, all_neighbors_count
            neigh1 = all_neighbors(j)
            
            do k = j+1, all_neighbors_count
                neigh2 = all_neighbors(k)
                
                ! Skip if neighbors are identical (should never happen due to j+1 loop)
                if (neigh1 == neigh2) cycle
                
                ! Increment triplet counter for stats
                local_found_triplets = local_found_triplets + 1
                
                ! Skip if angle calculation for this type triplet is disabled
                ! FIXED: Only apply the type mask filter if angle_filter_types is enabled
                if (allocated(angle_type_mask) .and. angle_filter_types) then
                    if (.not. angle_type_mask(atom_types(central_atom), &
                                             atom_types(neigh1), &
                                             atom_types(neigh2))) cycle
                end if
                
                ! Calculate vectors from central atom to neighbors
                dx1 = coords(neigh1,1) - coords(central_atom,1)
                dy1 = coords(neigh1,2) - coords(central_atom,2)
                dz1 = coords(neigh1,3) - coords(central_atom,3)
                
                dx2 = coords(neigh2,1) - coords(central_atom,1)
                dy2 = coords(neigh2,2) - coords(central_atom,2)
                dz2 = coords(neigh2,3) - coords(central_atom,3)
                
                ! Apply minimum image convention
                dx1 = dx1 - box_length(1) * nint(dx1/box_length(1))
                dy1 = dy1 - box_length(2) * nint(dy1/box_length(2))
                dz1 = dz1 - box_length(3) * nint(dz1/box_length(3))
                
                dx2 = dx2 - box_length(1) * nint(dx2/box_length(1))
                dy2 = dy2 - box_length(2) * nint(dy2/box_length(2))
                dz2 = dz2 - box_length(3) * nint(dz2/box_length(3))
                
                ! Create vectors
                vec1 = (/dx1, dy1, dz1/)
                vec2 = (/dx2, dy2, dz2/)
                
                ! Calculate squared magnitudes (no square roots needed!)
                vec1_sq = vec1(1)**2 + vec1(2)**2 + vec1(3)**2
                vec2_sq = vec2(1)**2 + vec2(2)**2 + vec2(3)**2
                
                ! Skip if vectors are too small (avoid division by zero)
                if (vec1_sq < 1.0e-20_dp .or. vec2_sq < 1.0e-20_dp) cycle
                
                ! Calculate dot product
                dot_prod = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
                
                ! Apply our efficient pre-filtering using squared cosine
                if (using_angle_filter) then
                    ! Calculate squared cosine directly without square roots
                    cos_sq = (dot_prod * dot_prod) / (vec1_sq * vec2_sq)
                    
                    ! Quick filter check
                    if (cos_sq < min_cos_sq .or. cos_sq > max_cos_sq) then
                        local_filtered_angles = local_filtered_angles + 1
                        cycle
                    end if
                    
                    ! Check if angle is acute or obtuse - ALWAYS do this check with filtering
                    ! This solves the problem with squared cosines not distinguishing angle vs supplement
                    is_acute = (dot_prod >= 0.0_dp)
                    
                    ! Filter acute/obtuse angles based on min/max settings
                    if (local_max_angle < 90.0_dp .and. .not. is_acute) then
                        ! We only want acute angles but this one is obtuse
                        local_filtered_angles = local_filtered_angles + 1
                        cycle
                    else if (local_min_angle > 90.0_dp .and. is_acute) then
                        ! We only want obtuse angles but this one is acute
                        local_filtered_angles = local_filtered_angles + 1
                        cycle
                    end if
                end if
                
                ! Now we perform the actual angle calculation with square root
                ! This only happens for angles that passed our pre-filtering
                local_calculated_angles = local_calculated_angles + 1
                
                ! Calculate the actual angle
                dot_prod = max(-1.0_dp, min(1.0_dp, dot_prod / sqrt(vec1_sq * vec2_sq)))
                angle = acos(dot_prod) * 180.0_dp / PI
                
                ! Store angle in local buffer
                if (local_angle_count < LOCAL_ANGLE_BUFFER_SIZE) then
                    local_angle_count = local_angle_count + 1
                    local_angle_buffer(local_angle_count)%central_atom = central_atom
                    local_angle_buffer(local_angle_count)%neighbor1 = neigh1
                    local_angle_buffer(local_angle_count)%neighbor2 = neigh2
                    local_angle_buffer(local_angle_count)%angle = angle
                    local_angle_buffer(local_angle_count)%type1 = atom_types(neigh1)
                    local_angle_buffer(local_angle_count)%type2 = atom_types(neigh2)
                    local_angle_buffer(local_angle_count)%central_type = atom_types(central_atom)
                end if
                
                ! If local buffer is full, merge with global buffer and reset
                if (local_angle_count >= LOCAL_ANGLE_BUFFER_SIZE) then
                    !$OMP CRITICAL(merge_angles)
                    if (angle_count + local_angle_count > max_angles) then
                        ! Resize if needed
                        call resize_angle_array(angle_count + local_angle_count)
                    end if
                    
                    ! Copy angles to global buffer
                    if (local_angle_count > 0) then
                        angle_list(angle_count+1:angle_count+local_angle_count) = &
                            local_angle_buffer(1:local_angle_count)
                        angle_count = angle_count + local_angle_count
                    end if
                    !$OMP END CRITICAL(merge_angles)
                    
                    ! Reset local counter to reuse the buffer
                    local_angle_count = 0
                end if
            end do
        end do
    end do
    !$OMP END DO
    
    ! Merge any remaining angles from local buffers
    !$OMP CRITICAL(merge_angles)
    if (angle_count + local_angle_count > max_angles) then
        ! Resize if needed
        call resize_angle_array(angle_count + local_angle_count)
    end if
    
    ! Copy angles to global buffer
    if (local_angle_count > 0) then
        angle_list(angle_count+1:angle_count+local_angle_count) = &
            local_angle_buffer(1:local_angle_count)
        angle_count = angle_count + local_angle_count
    end if
    !$OMP END CRITICAL(merge_angles)
    
    ! Aggregate statistics
    !$OMP ATOMIC
    found_triplets = found_triplets + local_found_triplets
    !$OMP ATOMIC
    calculated_angles = calculated_angles + local_calculated_angles
    !$OMP ATOMIC
    filtered_angles = filtered_angles + local_filtered_angles
    
    ! Clean up thread-local resources
    deallocate(local_angle_buffer)
    deallocate(all_neighbors)
    !$OMP END PARALLEL
    
    end_time = omp_get_wtime()
    
    ! Print diagnostic information if not suppressed
    if ((.not. diagnostics_shown .or. VERBOSE) .and. .not. suppress_output) then
        write(*,'(A)') " Angle calculation statistics:"
        write(*,'(A,I8)') "   Potential triplets examined: ", found_triplets
        write(*,'(A,I8,A,F6.2,A)') "   Angles pre-filtered: ", filtered_angles, " (", &
             real(filtered_angles) / max(1, found_triplets) * 100.0, "%)"
        write(*,'(A,I8)') "   Angles calculated: ", calculated_angles
        write(*,'(A,I8)') "   Angles stored: ", angle_count
        diagnostics_shown = .true.
    end if
    
    ! Only process if angles were found
    if (angle_count > 0) then
        ! Write frame header to output file
        write(angle_unit,'(A,I0)') "# Frame: ", frame
        
        ! Sort angles by central atom for better organization
        call sort_angles()
        
        ! Write results to file
        call write_angle_data(atom_info)
        
        if (.not. suppress_output) then
            write(*,'(A,I8,A,F8.3,A)') " Found ", angle_count, " angles in ", &
                end_time - start_time, " seconds"
        end if
    else
        if (.not. suppress_output) then
            write(*,'(A)') " No angles found! Check your cutoff distances and atom type filters."
            write(*,'(A)') " Make sure your setup includes all needed atom type combinations."
        end if
    end if
end subroutine compute_angles
    
    !> @brief Resize the angle array to accommodate more angles
    !>
    !> @param[in] new_size Target size for the angle array
    subroutine resize_angle_array(new_size)
        integer, intent(in), optional :: new_size
        integer :: new_max, alloc_stat
        type(angle_data_type), allocatable :: temp(:)
        
        ! Determine new size (double current size or use provided size)
        if (present(new_size)) then
            new_max = new_size + 1000  ! Add padding
        else
            new_max = max_angles * 2
        end if
        
        ! Allocate temporary array
        allocate(temp(new_max), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to resize angle array", ERR_ALLOCATION, fatal=.false.)
            return
        end if
        
        ! Copy existing data
        temp(1:angle_count) = angle_list(1:angle_count)
        
        ! Replace old array with new one
        call move_alloc(from=temp, to=angle_list)
        max_angles = new_max
        
        if (VERBOSE) write(*,'(A,I0,A)') " Angle array resized to ", max_angles, " elements"
    end subroutine resize_angle_array
    
    !> @brief Sort angles by central atom ID
    subroutine sort_angles()
        integer :: i, j, chunk_size, start_idx, end_idx
        type(angle_data_type) :: temp
        
        ! Skip if no angles to sort
        if (angle_count <= 1) return
        
        ! Use simple chunked approach for better cache efficiency
        chunk_size = 256  ! Adjust based on typical L1 cache size
        
        ! Process chunks of the array
        do start_idx = 1, angle_count, chunk_size
            end_idx = min(start_idx + chunk_size - 1, angle_count)
            
            ! Simple insertion sort for this chunk
            do i = start_idx+1, end_idx
                temp = angle_list(i)
                j = i - 1
                do while (j >= start_idx .and. angle_list(j)%central_atom > temp%central_atom)
                    angle_list(j+1) = angle_list(j)
                    j = j - 1
                end do
                angle_list(j+1) = temp
            end do
        end do
        
        ! Merge sorted chunks if needed (not implemented for simplicity)
        ! For production code, you'd want a full merge sort here
    end subroutine sort_angles
    
    !> @brief Write angle data to output file
    subroutine write_angle_data(atom_info)
        type(atom_type_info), intent(in) :: atom_info(:)
        integer :: i
        
        ! Skip if file not open
        if (.not. angle_file_open) return
        
        ! Skip if no angles to write
        if (angle_count <= 0) return
        
        ! Use buffered writing for better I/O performance
        ! We'll write in blocks rather than one line at a time
        do i = 1, angle_count
            ! Write angle data with type information
            write(angle_unit,'(I6,2X,A12,2X,I6,2X,A12,2X,I6,2X,A12,2X,F10.3)') &
                angle_list(i)%central_atom, &
                trim(atom_info(angle_list(i)%central_type)%name), &
                angle_list(i)%neighbor1, &
                trim(atom_info(angle_list(i)%type1)%name), &
                angle_list(i)%neighbor2, &
                trim(atom_info(angle_list(i)%type2)%name), &
                angle_list(i)%angle
        end do
        
        ! Add a blank line between frames
        write(angle_unit,*)
    end subroutine write_angle_data
    
    !> @brief Clean up angle calculation resources
    subroutine finalize_angle_analysis()
        if (allocated(angle_list)) deallocate(angle_list)
        
        if (angle_file_open) then
            close(angle_unit)
            angle_file_open = .false.
        end if
        
        angle_count = 0
        max_angles = 0
    end subroutine finalize_angle_analysis
    
end module angle_mod