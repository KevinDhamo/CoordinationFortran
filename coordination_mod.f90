!> @brief Module for calculating atomic coordination numbers
!>
!> This module implements the core algorithm for calculating atomic coordination
!> numbers between different atom types. It uses a cell-based spatial partitioning approach
!> for efficient neighbor finding, but can also use a direct calculation approach for 
!> validation purposes.
!>
!> The module supports parallel execution through OpenMP and can track both
!> coordination numbers and the specific neighbor atom IDs.
module coordination_mod
    use types_mod
    use config_mod
    use cell_list_mod
    use error_mod
    use atom_selection_mod
    use omp_lib
    implicit none
    private

    !===============================================================================
    ! Public procedures
    !===============================================================================
    public :: initialize_coordination
    public :: calculate_coordination
    public :: cleanup_coordination

    !===============================================================================
    ! Public variables
    !===============================================================================
    !> Coordination numbers matrix (atom_index, pair_index)
    public :: coord_numbers
    
    !> Atom type pair definitions
    public :: pairs
    
    !> Number of atom type pairs
    public :: n_pairs
    
    !> Atom neighbor tracking
    public :: coord_neighbors
    
    !> Flag to show neighbor atom IDs in output
    public :: show_neighbors

    !> Count of neighbors for each atom/pair
    public :: neighbor_counts
    
    !===============================================================================
    ! Module variables
    !===============================================================================
    !> Coordination numbers matrix (atom_index, pair_index)
    integer, allocatable :: coord_numbers(:,:)
    
    !> Atom type pair definitions
    type(pair_type), allocatable :: pairs(:)
    
    !> Number of atom type pairs
    integer :: n_pairs
    
    !===============================================================================
    ! Neighbor tracking variables
    !===============================================================================
    !> Flag to enable/disable neighbor tracking
    logical :: show_neighbors = .false.
    
    !> Maximum number of neighbors to track per atom/pair
    integer, parameter :: MAX_NEIGHBORS = 50
    
    !> Array of neighbor atom IDs (atom_index, pair_index, neighbor_index)
    integer, allocatable :: coord_neighbors(:,:,:)
    
    !> Count of neighbors for each atom/pair
    integer, allocatable :: neighbor_counts(:,:)
    
    !> Flag to track if setup information has been displayed
    logical :: setup_shown = .false.

contains
    !> @brief Initialize coordination analysis
    !>
    !> Sets up data structures for coordination analysis, including
    !> pair cutoff distances and coordination number arrays.
    !>
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] n_atom_types Number of atom types in the system
    !> @param[in] atom_info Array of atom type information
    subroutine initialize_coordination(n_atoms, n_atom_types, atom_info)
        integer, intent(in) :: n_atoms, n_atom_types
        type(atom_type_info), intent(in) :: atom_info(:)
        integer :: i, j, pair_idx, alloc_stat
        character(20) :: user_input
        real(dp) :: cutoff_value
        
        ! Calculate number of unique pairs
        n_pairs = 0
        do i = 1, n_atom_types
            if (.not. atom_info(i)%include) cycle
            do j = i, n_atom_types
                if (.not. atom_info(j)%include) cycle
                n_pairs = n_pairs + 1
            end do
        end do
        
        ! Allocate coordination number array
        allocate(coord_numbers(n_atoms, n_pairs), stat=alloc_stat)
        call check_allocation(alloc_stat, "coordination arrays")
        
        ! Allocate neighbor tracking arrays
        allocate(coord_neighbors(n_atoms, n_pairs, MAX_NEIGHBORS), &
                 neighbor_counts(n_atoms, n_pairs), stat=alloc_stat)
        call check_allocation(alloc_stat, "neighbor tracking arrays")
        
        ! Check if pairs is already allocated (from setup file)
        if (.not. allocated(pairs)) then
            allocate(pairs(n_pairs), stat=alloc_stat)
            call check_allocation(alloc_stat, "pairs array")
            
            ! Initialize pair information through user prompts
            write(*,'(A)') ' Setting up pair cutoff distances:'
            pair_idx = 1
            do i = 1, n_atom_types
                if (.not. atom_info(i)%include) cycle
                do j = i, n_atom_types
                    if (.not. atom_info(j)%include) cycle
                    
                    pairs(pair_idx)%type1_id = i
                    pairs(pair_idx)%type2_id = j
                    
                    write(*,'(A,I0,A,A,A,I0,A,A,A,F5.2,A)', advance='no') &
                        '   Enter cutoff for type ', i, ' (', trim(atom_info(i)%name), ') - type ', j, &
                        ' (', trim(atom_info(j)%name), ') (Å) [', DEFAULT_CUTOFF, ']: '
                    
                    read(*,'(A)') user_input
                    if (len_trim(user_input) == 0) then
                        cutoff_value = DEFAULT_CUTOFF
                    else
                        read(user_input,*,iostat=alloc_stat) cutoff_value
                        if (alloc_stat /= 0) then
                            call handle_error("Invalid cutoff value", ERR_INVALID_PARAM)
                        end if
                    end if
                    
                    ! Validate cutoff value
                    if (cutoff_value <= 0.0_dp) then
                        call handle_error("Cutoff value must be positive", ERR_INVALID_PARAM)
                    end if
                    
                    pairs(pair_idx)%cutoff = cutoff_value
                    pairs(pair_idx)%cutoff_sq = cutoff_value**2
                    pair_idx = pair_idx + 1
                end do
            end do
        end if
        
        ! Initialize coordination numbers and neighbor tracking
        coord_numbers = 0
        coord_neighbors = 0
        neighbor_counts = 0
        
        ! Check if neighbor tracking is enabled
        show_neighbors = show_atom_neighbors
    end subroutine initialize_coordination

    !> @brief Calculate coordination numbers for the current configuration
    !>
    !> Main function for calculating coordination numbers between all atoms
    !> in the system based on their positions and the specified cutoff distances.
    !> Can use either cell list or direct (O(n²)) calculation method based on 
    !> the use_cell_list configuration option.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
    !> @param[in] force_rebuild Optional flag to force a complete rebuild (default: .false.)
    subroutine calculate_coordination(coords, atom_types, n_atoms, box_length, frame, include_mask, force_rebuild)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        logical, intent(in), optional :: force_rebuild
        
        logical :: do_force_rebuild
        
        ! Process optional parameter
        do_force_rebuild = .false.
        if (present(force_rebuild)) do_force_rebuild = force_rebuild
        
        ! Use appropriate calculation method based on configuration
        if (use_cell_list) then
            ! Use cell list method for efficient O(n) calculation
            call calculate_coordination_cell(coords, atom_types, n_atoms, box_length, frame, include_mask, do_force_rebuild)
        else
            ! Use direct O(n²) method for validation
            call calculate_coordination_direct(coords, atom_types, n_atoms, box_length, include_mask)
        end if
    end subroutine calculate_coordination

    !> @brief Calculate coordination using cell lists
    !>
    !> Traditional coordination calculation using the cell list method.
    !> This processes atom pairs based on their cell assignments.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
    !> @param[in] force_rebuild Optional flag to force a complete rebuild (default: .false.)
    subroutine calculate_coordination_cell(coords, atom_types, n_atoms, box_length, frame, include_mask, force_rebuild)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        logical, intent(in), optional :: force_rebuild
        
        integer :: ix, iy, iz, nx, ny, nz
        integer :: neighbor_cells(27,3), n_neighbors
        logical :: do_force_rebuild
        
        ! Process optional parameter
        do_force_rebuild = .false.
        if (present(force_rebuild)) do_force_rebuild = force_rebuild
        
        ! Update cell lists - include ALL atoms of valid types (not just selected ones)
        call update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask, do_force_rebuild)
        
        ! Initialize coordination arrays
        coord_numbers = 0
        neighbor_counts = 0  ! Reset for this frame
        
        ! Get cell grid dimensions
        call get_cell_grid_dims(nx, ny, nz)
        
        ! Loop over all cells in parallel
        !$OMP PARALLEL DO PRIVATE(iy,iz,neighbor_cells,n_neighbors) SCHEDULE(dynamic)
        do ix = 1, nx
            do iy = 1, ny
                do iz = 1, nz
                    ! Get neighboring cells
                    call get_neighboring_cells(ix, iy, iz, neighbor_cells, n_neighbors)
                    
                    ! Process atoms in current cell
                    call process_cell_atoms(ix, iy, iz, coords, atom_types, &
                                          box_length, include_mask)
                                          
                    ! Process atoms in neighboring cells with ordering to avoid double counting
                    call process_neighbor_cells(neighbor_cells, n_neighbors, ix, iy, iz, &
                                              coords, atom_types, box_length, include_mask)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calculate_coordination_cell

    !> @brief Calculate coordination using direct O(n²) method
    !>
    !> This function calculates coordination numbers using a direct O(n²) approach
    !> without spatial partitioning. This is slower but useful for validation and testing.
    !> FIXED: Now correctly handles atom selection - only calculates coordination FOR 
    !> selected atoms, but considers ALL atoms as potential neighbors.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine calculate_coordination_direct(coords, atom_types, n_atoms, box_length, include_mask)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: i, j
        
        ! Reset coordination arrays
        coord_numbers = 0
        neighbor_counts = 0
        
        ! Display information about using direct method (only once)
        if (.not. setup_shown) then
            write(*,'(A)') " Using direct O(n²) calculation method for coordination"
            write(*,'(A)') " Note: This method is slower but can be used for validation"
            setup_shown = .true.
        end if
        
        ! FIXED: Only calculate coordination FOR selected atoms, but consider ALL atoms as neighbors
        !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(guided)
        do i = 1, n_atoms
            ! Skip if atom type is excluded OR atom is not selected
            if (.not. include_mask(atom_types(i)) .or. .not. is_atom_selected(i)) cycle
            
            ! Check this atom against ALL other atoms (not just selected ones)
            do j = 1, n_atoms
                ! Skip self-coordination
                if (i == j) cycle
                
                ! Skip if neighbor atom type is excluded (but don't check selection)
                if (.not. include_mask(atom_types(j))) cycle
                
                ! Process this atom pair - but only update coordination for selected atoms
                call process_atom_pair_selective(i, j, coords, atom_types, box_length, include_mask)
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calculate_coordination_direct

    !> @brief Process atom pairs within a single cell
    !>
    !> FIXED: Now properly handles atom selection in cell-based method.
    !>
    !> @param[in] ix X index of the cell
    !> @param[in] iy Y index of the cell
    !> @param[in] iz Z index of the cell
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine process_cell_atoms(ix, iy, iz, coords, atom_types, box_length, include_mask)
        integer, intent(in) :: ix, iy, iz
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: i, j, atom_i, atom_j
        type(cell_type) :: current_cell
        
        ! Get current cell
        call get_cell(ix, iy, iz, current_cell)
        
        ! Process atom pairs within the cell
        do i = 1, current_cell%n_atoms
            atom_i = current_cell%atoms(i)
            ! Skip if atom type is excluded OR atom is not selected
            if (.not. include_mask(atom_types(atom_i)) .or. .not. is_atom_selected(atom_i)) cycle
            
            do j = 1, current_cell%n_atoms
                ! Skip self-coordination
                if (i == j) cycle
                
                atom_j = current_cell%atoms(j)
                ! Only check if neighbor atom type is included (don't check selection)
                if (.not. include_mask(atom_types(atom_j))) cycle
                
                call process_atom_pair_selective(atom_i, atom_j, coords, atom_types, &
                                               box_length, include_mask)
            end do
        end do
    end subroutine process_cell_atoms

    !> @brief Process atom pairs between the current cell and neighboring cells
    !>
    !> FIXED: Now properly handles atom selection in cell-based method.
    !>
    !> @param[in] neighbor_cells Array of neighbor cell indices
    !> @param[in] n_neighbors Number of neighboring cells
    !> @param[in] ix X index of the current cell
    !> @param[in] iy Y index of the current cell
    !> @param[in] iz Z index of the current cell
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine process_neighbor_cells(neighbor_cells, n_neighbors, ix, iy, iz, &
                                    coords, atom_types, box_length, include_mask)
        integer, intent(in) :: neighbor_cells(:,:), n_neighbors, ix, iy, iz
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: k, current_cell_idx, neighbor_cell_idx
        
        ! Calculate cell index for the current cell
        current_cell_idx = cell_index(ix, iy, iz)
        
        do k = 1, n_neighbors
            ! Calculate cell index for the neighbor cell
            neighbor_cell_idx = cell_index(neighbor_cells(k,1), neighbor_cells(k,2), neighbor_cells(k,3))
            
            ! Only process if current cell comes before neighbor cell
            ! This prevents double counting pairs across different cells
            if (current_cell_idx < neighbor_cell_idx) then
                call process_cell_pair(ix, iy, iz, &
                                    neighbor_cells(k,1), neighbor_cells(k,2), neighbor_cells(k,3), &
                                    coords, atom_types, box_length, include_mask)
            end if
        end do
    end subroutine process_neighbor_cells

    !> @brief Process atoms between two different cells
    !>
    !> FIXED: Now properly handles atom selection in cell-based method.
    !>
    !> @param[in] ix1 X index of first cell
    !> @param[in] iy1 Y index of first cell
    !> @param[in] iz1 Z index of first cell
    !> @param[in] ix2 X index of second cell
    !> @param[in] iy2 Y index of second cell
    !> @param[in] iz2 Z index of second cell
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine process_cell_pair(ix1, iy1, iz1, ix2, iy2, iz2, &
                                coords, atom_types, box_length, include_mask)
        integer, intent(in) :: ix1, iy1, iz1, ix2, iy2, iz2
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        integer :: i, j, atom_i, atom_j
        type(cell_type) :: cell1, cell2
        
        ! Skip if same cell (already processed in process_cell_atoms)
        if (ix1 == ix2 .and. iy1 == iy2 .and. iz1 == iz2) return
        
        ! Get cells
        call get_cell(ix1, iy1, iz1, cell1)
        call get_cell(ix2, iy2, iz2, cell2)
        
        ! Process atoms from cell1 to cell2
        do i = 1, cell1%n_atoms
            atom_i = cell1%atoms(i)
            ! Only process if atom is selected and type is included
            if (.not. include_mask(atom_types(atom_i)) .or. .not. is_atom_selected(atom_i)) cycle
            
            do j = 1, cell2%n_atoms
                atom_j = cell2%atoms(j)
                ! Only check if neighbor atom type is included (don't check selection)
                if (.not. include_mask(atom_types(atom_j))) cycle
                
                call process_atom_pair_selective(atom_i, atom_j, coords, atom_types, &
                                               box_length, include_mask)
            end do
        end do
        
        ! Process atoms from cell2 to cell1 (reverse direction)
        do j = 1, cell2%n_atoms
            atom_j = cell2%atoms(j)
            ! Only process if atom is selected and type is included
            if (.not. include_mask(atom_types(atom_j)) .or. .not. is_atom_selected(atom_j)) cycle
            
            do i = 1, cell1%n_atoms
                atom_i = cell1%atoms(i)
                ! Only check if neighbor atom type is included (don't check selection)
                if (.not. include_mask(atom_types(atom_i))) cycle
                
                call process_atom_pair_selective(atom_j, atom_i, coords, atom_types, &
                                               box_length, include_mask)
            end do
        end do
    end subroutine process_cell_pair

    !> @brief Process a single atom pair to check for coordination (legacy version)
    !>
    !> DEPRECATED: This is the old version that incorrectly filtered both atoms.
    !> Use process_atom_pair_selective instead.
    !>
    !> @param[in] atom_i Index of first atom
    !> @param[in] atom_j Index of second atom
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine process_atom_pair(atom_i, atom_j, coords, atom_types, box_length, include_mask)
        integer, intent(in) :: atom_i, atom_j
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        ! Redirect to the new selective version
        call process_atom_pair_selective(atom_i, atom_j, coords, atom_types, box_length, include_mask)
    end subroutine process_atom_pair

    !> @brief Process a single atom pair with proper atom selection handling
    !>
    !> FIXED: Core function that calculates if two atoms are coordinated based on
    !> their distance and the cutoff for their atom types. Only updates coordination
    !> numbers for SELECTED atoms, but considers all atoms as potential neighbors.
    !>
    !> @param[in] atom_i Index of first atom
    !> @param[in] atom_j Index of second atom
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine process_atom_pair_selective(atom_i, atom_j, coords, atom_types, box_length, include_mask)
        integer, intent(in) :: atom_i, atom_j
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        
        real(dp) :: dx, dy, dz, r_sq
        integer :: pair_idx
        logical :: atom_i_selected, atom_j_selected
        
        ! Skip if either atom type is excluded from analysis
        if (.not. include_mask(atom_types(atom_i)) .or. .not. include_mask(atom_types(atom_j))) return
        
        ! Check which atoms are selected
        atom_i_selected = is_atom_selected(atom_i)
        atom_j_selected = is_atom_selected(atom_j)
        
        ! Skip if neither atom is selected (no point in calculating)
        if (.not. atom_i_selected .and. .not. atom_j_selected) return
        
        ! Calculate distance vector between atoms
        dx = coords(atom_i,1) - coords(atom_j,1)
        dy = coords(atom_i,2) - coords(atom_j,2)
        dz = coords(atom_i,3) - coords(atom_j,3)
        
        ! Apply minimum image convention for periodic boundaries
        dx = dx - nint(dx/box_length(1)) * box_length(1)
        dy = dy - nint(dy/box_length(2)) * box_length(2)
        dz = dz - nint(dz/box_length(3)) * box_length(3)
        
        ! Calculate squared distance
        r_sq = dx*dx + dy*dy + dz*dz
        
        ! Find pair index and update coordination if within cutoff
        call get_pair_index(atom_types(atom_i), atom_types(atom_j), pair_idx)
        if (pair_idx > 0 .and. r_sq <= pairs(pair_idx)%cutoff_sq) then
            ! Update coordination numbers ONLY for selected atoms
            if (atom_i_selected) then
                !$OMP ATOMIC
                coord_numbers(atom_i,pair_idx) = coord_numbers(atom_i,pair_idx) + 1
            end if
            
            if (atom_j_selected) then
                !$OMP ATOMIC
                coord_numbers(atom_j,pair_idx) = coord_numbers(atom_j,pair_idx) + 1
            end if
            
            ! Update neighbor lists if enabled (only for selected atoms)
            if (show_neighbors) then
                call update_neighbor_lists_selective(atom_i, atom_j, pair_idx, atom_i_selected, atom_j_selected)
            end if
        end if
    end subroutine process_atom_pair_selective

    !> @brief Update neighbor tracking lists for coordination visualization (legacy version)
    !>
    !> DEPRECATED: Use update_neighbor_lists_selective instead.
    !>
    !> @param[in] atom_i Index of first atom
    !> @param[in] atom_j Index of second atom
    !> @param[in] pair_idx Index of the atom type pair
    subroutine update_neighbor_lists(atom_i, atom_j, pair_idx)
        integer, intent(in) :: atom_i, atom_j, pair_idx
        
        ! Redirect to selective version
        call update_neighbor_lists_selective(atom_i, atom_j, pair_idx, .true., .true.)
    end subroutine update_neighbor_lists

    !> @brief Update neighbor tracking lists with proper atom selection handling
    !>
    !> FIXED: Updates the neighbor tracking arrays used for visualizing coordination.
    !> Only updates neighbor lists for atoms that are actually selected.
    !>
    !> @param[in] atom_i Index of first atom
    !> @param[in] atom_j Index of second atom
    !> @param[in] pair_idx Index of the atom type pair
    !> @param[in] atom_i_selected Whether atom_i is selected
    !> @param[in] atom_j_selected Whether atom_j is selected
    subroutine update_neighbor_lists_selective(atom_i, atom_j, pair_idx, atom_i_selected, atom_j_selected)
        integer, intent(in) :: atom_i, atom_j, pair_idx
        logical, intent(in) :: atom_i_selected, atom_j_selected
        integer :: count_i, count_j
        
        ! Update neighbor counts and lists only for selected atoms
        if (atom_i_selected) then
            !$OMP ATOMIC
            neighbor_counts(atom_i,pair_idx) = neighbor_counts(atom_i,pair_idx) + 1
            
            !$OMP CRITICAL(get_count_i)
            count_i = neighbor_counts(atom_i,pair_idx)
            !$OMP END CRITICAL(get_count_i)
            
            ! Store neighbor if we haven't exceeded the maximum
            if (count_i <= MAX_NEIGHBORS) then
                !$OMP CRITICAL(update_neighbors_i)
                coord_neighbors(atom_i,pair_idx,count_i) = atom_j
                !$OMP END CRITICAL(update_neighbors_i)
            end if
        end if
        
        if (atom_j_selected) then
            !$OMP ATOMIC
            neighbor_counts(atom_j,pair_idx) = neighbor_counts(atom_j,pair_idx) + 1
            
            !$OMP CRITICAL(get_count_j)
            count_j = neighbor_counts(atom_j,pair_idx)
            !$OMP END CRITICAL(get_count_j)
            
            ! Store neighbor if we haven't exceeded the maximum
            if (count_j <= MAX_NEIGHBORS) then
                !$OMP CRITICAL(update_neighbors_j)
                coord_neighbors(atom_j,pair_idx,count_j) = atom_i
                !$OMP END CRITICAL(update_neighbors_j)
            end if
        end if
    end subroutine update_neighbor_lists_selective

    !> @brief Get the pair index for two atom types
    !>
    !> Returns the index in the pairs array for the given atom types,
    !> taking care to handle the ordering (always using the lower type ID first).
    !>
    !> @param[in] type1 Type ID of first atom
    !> @param[in] type2 Type ID of second atom
    !> @param[out] idx Index of the pair in the pairs array (-1 if not found)
    subroutine get_pair_index(type1, type2, idx)
        integer, intent(in) :: type1, type2
        integer, intent(out) :: idx
        integer :: t1, t2
        
        ! Order types for consistent pair lookup (smaller ID first)
        t1 = min(type1, type2)
        t2 = max(type1, type2)
        
        ! Find matching pair in pairs array
        do idx = 1, n_pairs
            if (pairs(idx)%type1_id == t1 .and. pairs(idx)%type2_id == t2) return
        end do
        
        ! Pair not found
        idx = -1
    end subroutine get_pair_index

    !> @brief Helper function to calculate a unique cell index
    !>
    !> Generates a unique integer identifier for a cell based on its 3D indices
    !> This is used for consistent ordering of cell pairs to avoid double counting
    !>
    !> @param[in] ix X-coordinate of the cell
    !> @param[in] iy Y-coordinate of the cell
    !> @param[in] iz Z-coordinate of the cell
    !> @return Unique integer identifier for the cell
    function cell_index(ix, iy, iz) result(idx)
        integer, intent(in) :: ix, iy, iz
        integer :: idx
        integer :: nx, ny, nz
        
        ! Get grid dimensions
        call get_cell_grid_dims(nx, ny, nz)
        
        ! Calculate unique index based on cell coordinates
        idx = ix + (iy-1)*nx + (iz-1)*nx*ny
    end function cell_index

    !> @brief Clean up coordination analysis resources
    !>
    !> Deallocates all arrays used by the coordination module.
    subroutine cleanup_coordination()
        if (allocated(coord_numbers)) deallocate(coord_numbers)
        if (allocated(coord_neighbors)) deallocate(coord_neighbors)
        if (allocated(neighbor_counts)) deallocate(neighbor_counts)
        if (allocated(pairs)) deallocate(pairs)
        
        setup_shown = .false.  ! Reset for potential reuse
    end subroutine cleanup_coordination

end module coordination_mod