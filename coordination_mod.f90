!> @brief Module for calculating atomic coordination numbers
!>
!> This module implements the core algorithm for calculating atomic coordination
!> numbers between different atom types. It supports two methods for neighbor finding:
!> 1. Cell List Method (traditional): Uses a cell-based spatial partitioning approach
!> 2. Verlet List Method (optional): Uses pre-computed neighbor lists with a skin distance
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
    
    !===============================================================================
    ! Verlet list variables
    !===============================================================================
    !> Array of Verlet lists for each atom
    type(verlet_list_type), allocatable :: verlet_lists(:)
    
    !> Previous atom coordinates when Verlet lists were built
    real(dp), allocatable :: verlet_reference_coords(:,:)
    
    !> Maximum displacement before Verlet lists need rebuilding
    real(dp) :: verlet_max_displacement
    
    !> Frame number when Verlet lists were last built
    integer :: verlet_last_build_frame = -1
    
    !> Total distance cutoff for Verlet lists (cutoff + skin)
    real(dp) :: verlet_cutoff_distance
    
    !> Flag indicating if Verlet lists are in use
    logical :: verlet_lists_active = .false.
    
    !> Flag indicating if Verlet lists need rebuilding
    logical :: verlet_lists_need_rebuild = .true.
    
    !> Flag to track if Verlet info has been displayed
    logical :: verlet_info_shown = .false.
    
    !> Type for Verlet list statistics reporting
    type :: verlet_stats_type
        integer :: total_builds = 0          !< Total number of Verlet list builds
        integer :: frames_processed = 0      !< Total frames processed
        real(dp) :: build_time = 0.0_dp      !< Total time spent building lists (seconds)
        real(dp) :: avg_neighbors = 0.0_dp   !< Average number of neighbors per atom
        integer :: min_neighbors = 0         !< Minimum neighbors for any atom
        integer :: max_neighbors = 0         !< Maximum neighbors for any atom
        real(dp) :: avg_build_time = 0.0_dp  !< Average time per build (seconds)
    end type verlet_stats_type
    
    !> Verlet list statistics
    type(verlet_stats_type) :: verlet_stats
    
    !> Make verlet_stats available to output module
    public :: verlet_stats

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
        
        ! Initialize Verlet list statistics
        verlet_stats = verlet_stats_type()
        
        ! Initialize Verlet lists if enabled
        if (use_verlet_lists) then
            call initialize_verlet_lists(n_atoms)
        end if
    end subroutine initialize_coordination
    
    !> @brief Initialize Verlet lists for all atoms
    !>
    !> Sets up Verlet list data structures if Verlet list mode is enabled.
    !> Calculates the maximum cutoff plus skin distance and allocates memory.
    !>
    !> @param[in] n_atoms Number of atoms in the system
    subroutine initialize_verlet_lists(n_atoms)
        integer, intent(in) :: n_atoms
        integer :: i, alloc_stat
        real(dp) :: max_cutoff
        
        ! Find maximum cutoff distance among all pairs
        max_cutoff = 0.0_dp
        do i = 1, n_pairs
            max_cutoff = max(max_cutoff, pairs(i)%cutoff)
        end do
        
        ! Calculate Verlet cutoff distance (cutoff + skin)
        verlet_cutoff_distance = max_cutoff + verlet_skin_distance
        
        ! Calculate maximum displacement before rebuild needed
        ! Typically half the skin distance is used
        verlet_max_displacement = verlet_skin_distance * verlet_rebuild_threshold
        
        ! Allocate Verlet lists
        if (allocated(verlet_lists)) deallocate(verlet_lists)
        allocate(verlet_lists(n_atoms), stat=alloc_stat)
        call check_allocation(alloc_stat, "verlet lists")
        
        ! Initialize each Verlet list
        do i = 1, n_atoms
            call verlet_lists(i)%init(INITIAL_VERLET_CAPACITY)
        end do
        
        ! Allocate reference coordinates array
        if (allocated(verlet_reference_coords)) deallocate(verlet_reference_coords)
        allocate(verlet_reference_coords(n_atoms, 3), stat=alloc_stat)
        call check_allocation(alloc_stat, "verlet reference coordinates")
        
        ! Mark Verlet lists as active but needing rebuild
        verlet_lists_active = .true.
        verlet_lists_need_rebuild = .true.
        verlet_last_build_frame = -1
        
        ! Display information if not already shown
        if (.not. verlet_info_shown .and. VERBOSE) then
            write(*,'(A)') ' Verlet list mode initialized:'
            write(*,'(A,F10.4,A)') '   Cutoff + skin distance: ', verlet_cutoff_distance, ' Å'
            write(*,'(A,F10.4,A)') '   Skin distance: ', verlet_skin_distance, ' Å'
            write(*,'(A,F10.4,A)') '   Rebuild threshold: ', verlet_max_displacement, ' Å'
            write(*,'(A,I0)') '   Initial neighbor capacity: ', INITIAL_VERLET_CAPACITY
            verlet_info_shown = .true.  ! Mark that we've shown the info
        end if
    end subroutine initialize_verlet_lists

    !> @brief Calculate coordination numbers for the current configuration
    !>
    !> Main function for calculating coordination numbers between all atoms
    !> in the system based on their positions and the specified cutoff distances.
    !> Uses either cell list or Verlet list approach based on configuration.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine calculate_coordination(coords, atom_types, n_atoms, box_length, frame, include_mask)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        
        ! If Verlet lists are active, use them for coordination calculation
        if (verlet_lists_active) then
            call calculate_coordination_verlet(coords, atom_types, n_atoms, box_length, frame, include_mask)
        else
            ! Otherwise use traditional cell list method
            call calculate_coordination_cell(coords, atom_types, n_atoms, box_length, frame, include_mask)
        end if
    end subroutine calculate_coordination
    
    !> @brief Calculate coordination using Verlet lists
    !>
    !> Calculates coordination numbers using pre-computed Verlet lists,
    !> updating the lists if necessary based on atom displacements.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
subroutine calculate_coordination_verlet(coords, atom_types, n_atoms, box_length, frame, include_mask)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: box_length(3)
    integer, intent(in) :: frame
    logical, intent(in) :: include_mask(:)
    
    integer :: i, j, neighbor_idx, atom_j, pair_idx
    real(dp) :: max_displacement, start_time, end_time
    real(dp) :: dx, dy, dz, r_sq
    
    ! Increment frame counter in statistics
    verlet_stats%frames_processed = verlet_stats%frames_processed + 1
    
    ! Check if we need to rebuild Verlet lists
    if (verlet_lists_need_rebuild .or. verlet_last_build_frame < 0) then
        ! Forced rebuild
        call build_verlet_lists(coords, atom_types, n_atoms, box_length, frame, include_mask)
    else
        ! Calculate maximum atom displacement since last Verlet list build
        max_displacement = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(dx,dy,dz) REDUCTION(max:max_displacement) SCHEDULE(guided)
        do i = 1, n_atoms
            ! Skip atoms of excluded types
            if (.not. include_mask(atom_types(i))) cycle
            
            ! Calculate displacement vector
            dx = coords(i,1) - verlet_reference_coords(i,1)
            dy = coords(i,2) - verlet_reference_coords(i,2)
            dz = coords(i,3) - verlet_reference_coords(i,3)
            
            ! Apply minimum image convention for periodic boundaries
            dx = dx - nint(dx/box_length(1)) * box_length(1)
            dy = dy - nint(dy/box_length(2)) * box_length(2)
            dz = dz - nint(dz/box_length(3)) * box_length(3)
            
            ! Calculate displacement magnitude
            max_displacement = max(max_displacement, sqrt(dx*dx + dy*dy + dz*dz))
        end do
        !$OMP END PARALLEL DO
        
        ! Check if displacement exceeds rebuild threshold
        if (max_displacement > verlet_max_displacement) then
            call build_verlet_lists(coords, atom_types, n_atoms, box_length, frame, include_mask)
        end if
    end if
    
    ! Initialize coordination number arrays
    coord_numbers = 0
    neighbor_counts = 0  ! Reset for this frame
    
    ! Process atom pairs from Verlet lists to calculate coordination numbers
    !$OMP PARALLEL DO PRIVATE(j,neighbor_idx,atom_j,pair_idx,dx,dy,dz,r_sq) SCHEDULE(guided)
    do i = 1, n_atoms
        ! Skip atoms that are not selected or excluded types
        if (.not. include_mask(atom_types(i)) .or. .not. is_atom_selected(i)) cycle
        
        ! Process each neighbor in the Verlet list
        do neighbor_idx = 1, verlet_lists(i)%n_neighbors
            atom_j = verlet_lists(i)%neighbors(neighbor_idx)
            pair_idx = verlet_lists(i)%pair_types(neighbor_idx)
            
            ! Skip invalid pairs or unselected atoms
            if (pair_idx <= 0 .or. .not. is_atom_selected(atom_j)) cycle
            
            ! Calculate current distance to check against actual cutoff
            dx = coords(i,1) - coords(atom_j,1)
            dy = coords(i,2) - coords(atom_j,2)
            dz = coords(i,3) - coords(atom_j,3)
            
            ! Apply minimum image convention for periodic boundaries
            dx = dx - nint(dx/box_length(1)) * box_length(1)
            dy = dy - nint(dy/box_length(2)) * box_length(2)
            dz = dz - nint(dz/box_length(3)) * box_length(3)
            
            ! Calculate squared distance
            r_sq = dx*dx + dy*dy + dz*dz
            
            ! Check if within actual cutoff (not skin distance)
            if (r_sq <= pairs(pair_idx)%cutoff_sq) then
                ! To avoid double counting, increment only when i < atom_j
                if (i < atom_j) then
                    ! Increment coordination numbers for both atoms
                    !$OMP ATOMIC
                    coord_numbers(i, pair_idx) = coord_numbers(i, pair_idx) + 1
                    !$OMP ATOMIC
                    coord_numbers(atom_j, pair_idx) = coord_numbers(atom_j, pair_idx) + 1
                    
                    ! Update neighbor lists if enabled
                    if (show_neighbors) then
                        call update_neighbor_lists(i, atom_j, pair_idx)
                    end if
                end if
            end if
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine calculate_coordination_verlet
    !> @brief Build Verlet lists for all atoms
    !>
    !> Constructs Verlet lists for all atoms using the cell list method
    !> to efficiently find atoms within the cutoff plus skin distance.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
subroutine build_verlet_lists(coords, atom_types, n_atoms, box_length, frame, include_mask)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: box_length(3)
    integer, intent(in) :: frame
    logical, intent(in) :: include_mask(:)
    
    integer :: i, max_len
    integer, allocatable :: temp_neighbors(:), temp_pair_types(:)
    integer :: n_neighbors, alloc_stat
    real(dp) :: start_time, end_time
    real(dp) :: total_neighbors
    integer :: local_min_neighbors, local_max_neighbors
    
    ! Start timing the build process
    start_time = omp_get_wtime()
    
    ! Allocate temporary arrays for neighbor finding
    max_len = 1000  ! Initial size, will grow if needed
    allocate(temp_neighbors(max_len), temp_pair_types(max_len), stat=alloc_stat)
    call check_allocation(alloc_stat, "temp neighbor arrays")
    
    ! Update cell list to ensure it's current
    call update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask)
    
    ! Save reference coordinates for displacement tracking
    verlet_reference_coords = coords
    
    ! Total neighbors for calculating average
    total_neighbors = 0.0_dp
    
    ! Initialize min/max neighbors
    local_min_neighbors = huge(local_min_neighbors)
    local_max_neighbors = 0
    
    ! Build Verlet lists for all atoms
    !$OMP PARALLEL DO PRIVATE(i,n_neighbors,temp_neighbors,temp_pair_types,max_len,alloc_stat) &
    !$OMP& REDUCTION(+:total_neighbors) &
    !$OMP& REDUCTION(min:local_min_neighbors) &
    !$OMP& REDUCTION(max:local_max_neighbors) &
    !$OMP& SCHEDULE(guided)
    do i = 1, n_atoms
        ! Skip atoms of excluded types
        if (.not. include_mask(atom_types(i))) cycle
        
        ! Clear Verlet list for this atom
        call verlet_lists(i)%cleanup()
        call verlet_lists(i)%init(INITIAL_VERLET_CAPACITY)
        
        ! Find all neighbors within cutoff+skin distance
        call find_neighbors_within_cutoff(i, coords, atom_types, box_length, &
                                      verlet_cutoff_distance, temp_neighbors, temp_pair_types, &
                                      n_neighbors, max_len, pairs, n_pairs)
        
        ! Store in Verlet list
        if (n_neighbors > 0) then
            ! Ensure Verlet list has enough capacity
            if (n_neighbors > size(verlet_lists(i)%neighbors)) then
                call verlet_lists(i)%resize(n_neighbors)
            end if
            
            ! Copy neighbors
            verlet_lists(i)%neighbors(1:n_neighbors) = temp_neighbors(1:n_neighbors)
            verlet_lists(i)%pair_types(1:n_neighbors) = temp_pair_types(1:n_neighbors)
            verlet_lists(i)%n_neighbors = n_neighbors
            
            ! Update statistics
            total_neighbors = total_neighbors + n_neighbors
            local_min_neighbors = min(local_min_neighbors, n_neighbors)
            local_max_neighbors = max(local_max_neighbors, n_neighbors)
        else
            verlet_lists(i)%n_neighbors = 0
        end if
    end do
    !$OMP END PARALLEL DO
    
    ! Clean up temporary arrays
    if (allocated(temp_neighbors)) deallocate(temp_neighbors)
    if (allocated(temp_pair_types)) deallocate(temp_pair_types)
    
    ! Update statistics
    if (n_atoms > 0) then
        verlet_stats%avg_neighbors = total_neighbors / real(n_atoms, dp)
    else
        verlet_stats%avg_neighbors = 0.0_dp
    end if
    
    ! Correct if no atoms had neighbors
    if (local_min_neighbors == huge(local_min_neighbors)) then
        local_min_neighbors = 0
    end if
    
    ! Update the module-level min/max statistics
    verlet_stats%min_neighbors = local_min_neighbors
    verlet_stats%max_neighbors = local_max_neighbors
    
    ! Record build time and update statistics
    end_time = omp_get_wtime()
    verlet_stats%build_time = verlet_stats%build_time + (end_time - start_time)
    verlet_stats%total_builds = verlet_stats%total_builds + 1
    
    if (verlet_stats%total_builds > 0) then
        verlet_stats%avg_build_time = verlet_stats%build_time / real(verlet_stats%total_builds, dp)
    end if
    
    ! Mark Verlet lists as valid and record frame
    verlet_lists_need_rebuild = .false.
    verlet_last_build_frame = frame
end subroutine build_verlet_lists
    
    !> @brief Update neighbor tracking lists for coordination visualization
    !>
    !> Updates the neighbor tracking arrays used for visualizing coordination.
    !> This is only used when show_neighbors is enabled and is separate from
    !> the Verlet lists used for calculation.
    !>
    !> @param[in] atom_i Index of first atom
    !> @param[in] atom_j Index of second atom
    !> @param[in] pair_idx Index of the atom type pair
    subroutine update_neighbor_lists(atom_i, atom_j, pair_idx)
        integer, intent(in) :: atom_i, atom_j, pair_idx
        integer :: count_i, count_j
        
        ! Update neighbor counts atomically
        !$OMP ATOMIC
        neighbor_counts(atom_i,pair_idx) = neighbor_counts(atom_i,pair_idx) + 1
        
        !$OMP ATOMIC
        neighbor_counts(atom_j,pair_idx) = neighbor_counts(atom_j,pair_idx) + 1
        
        ! Get the current counts after the atomic updates
        !$OMP CRITICAL(get_counts)
        count_i = neighbor_counts(atom_i,pair_idx)
        count_j = neighbor_counts(atom_j,pair_idx)
        !$OMP END CRITICAL(get_counts)
        
        ! Store neighbors if we haven't exceeded the maximum
        if (count_i <= MAX_NEIGHBORS) then
            !$OMP CRITICAL(update_neighbors_i)
            coord_neighbors(atom_i,pair_idx,count_i) = atom_j
            !$OMP END CRITICAL(update_neighbors_i)
        end if
        
        if (count_j <= MAX_NEIGHBORS) then
            !$OMP CRITICAL(update_neighbors_j)
            coord_neighbors(atom_j,pair_idx,count_j) = atom_i
            !$OMP END CRITICAL(update_neighbors_j)
        end if
    end subroutine update_neighbor_lists

    !> @brief Calculate coordination using cell lists
    !>
    !> Traditional coordination calculation using the cell list method.
    !> This is the original implementation that processes atom pairs
    !> based on their cell assignments.
    !>
    !> @param[in] coords Atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
    subroutine calculate_coordination_cell(coords, atom_types, n_atoms, box_length, frame, include_mask)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        
        integer :: ix, iy, iz, nx, ny, nz
        integer :: neighbor_cells(27,3), n_neighbors
        
        ! Update cell lists
        call update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask)
        
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
                                          
                    ! Process atoms in neighboring cells
                    call process_neighbor_cells(neighbor_cells, n_neighbors, ix, iy, iz, &
                                              coords, atom_types, box_length, include_mask)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calculate_coordination_cell

    !> @brief Process atom pairs within a single cell
    !>
    !> Computes coordination numbers for atoms located in the same cell.
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
            ! Check if atom is selected and included
            if (.not. include_mask(atom_types(atom_i)) .or. .not. is_atom_selected(atom_i)) cycle
            
            do j = i+1, current_cell%n_atoms
                atom_j = current_cell%atoms(j)
                ! Second atom is checked in process_atom_pair
                call process_atom_pair(atom_i, atom_j, coords, atom_types, &
                                     box_length, include_mask)
            end do
        end do
    end subroutine process_cell_atoms

    !> @brief Process atom pairs between the current cell and neighboring cells
    !>
    !> Processes coordination between atoms in the specified cell and its neighbors.
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
        
        integer :: k
        
        do k = 1, n_neighbors
            call process_cell_pair(ix, iy, iz, &
                                 neighbor_cells(k,1), neighbor_cells(k,2), neighbor_cells(k,3), &
                                 coords, atom_types, box_length, include_mask)
        end do
    end subroutine process_neighbor_cells

    !> @brief Process atoms between two different cells
    !>
    !> Calculates coordination numbers between atoms in two different cells.
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
        
        ! Process atoms between cells
        do i = 1, cell1%n_atoms
            atom_i = cell1%atoms(i)
            ! Check if atom is selected and included
            if (.not. include_mask(atom_types(atom_i)) .or. .not. is_atom_selected(atom_i)) cycle
            
            do j = 1, cell2%n_atoms
                atom_j = cell2%atoms(j)
                ! Second atom is checked in process_atom_pair
                call process_atom_pair(atom_i, atom_j, coords, atom_types, &
                                     box_length, include_mask)
            end do
        end do
    end subroutine process_cell_pair

    !> @brief Process a single atom pair to check for coordination
    !>
    !> Core function that calculates if two atoms are coordinated based on
    !> their distance and the cutoff for their atom types. Updates coordination
    !> numbers and neighbor lists.
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
        
        real(dp) :: dx, dy, dz, r_sq
        integer :: pair_idx
        
        ! Skip if second atom is not selected or type is excluded
        if (.not. include_mask(atom_types(atom_j)) .or. .not. is_atom_selected(atom_j)) return
        
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
            ! Update coordination numbers - use atomic operations for thread safety
            !$OMP ATOMIC
            coord_numbers(atom_i,pair_idx) = coord_numbers(atom_i,pair_idx) + 1
            !$OMP ATOMIC
            coord_numbers(atom_j,pair_idx) = coord_numbers(atom_j,pair_idx) + 1
            
            ! Update neighbor lists if enabled
            if (show_neighbors) then
                call update_neighbor_lists(atom_i, atom_j, pair_idx)
            end if
        end if
    end subroutine process_atom_pair

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

    !> @brief Clean up coordination analysis resources
    !>
    !> Deallocates all arrays used by the coordination module.
    subroutine cleanup_coordination()
        integer :: i
        
        if (allocated(coord_numbers)) deallocate(coord_numbers)
        if (allocated(coord_neighbors)) deallocate(coord_neighbors)
        if (allocated(neighbor_counts)) deallocate(neighbor_counts)
        if (allocated(pairs)) deallocate(pairs)
        
        ! Clean up Verlet lists if active
        if (verlet_lists_active) then
            if (allocated(verlet_lists)) then
                do i = 1, size(verlet_lists)
                    call verlet_lists(i)%cleanup()
                end do
                deallocate(verlet_lists)
            end if
            
            if (allocated(verlet_reference_coords)) deallocate(verlet_reference_coords)
            
            verlet_lists_active = .false.
        end if
        
        setup_shown = .false.  ! Reset for potential reuse
        verlet_info_shown = .false.
    end subroutine cleanup_coordination

end module coordination_mod