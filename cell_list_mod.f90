!> @brief Module implementing cell list algorithm for spatial partitioning
!>
!> This module provides a cell-based spatial partitioning system for efficient 
!> neighbor finding in molecular simulations. It divides the simulation box into
!> a grid of cells and assigns atoms to these cells, allowing for O(1) lookup of
!> neighboring atoms instead of O(N²) comparisons. The module supports dynamic
!> updates, selective rebuilding, parallel execution with OpenMP, and helper
!> functions for Verlet list construction.
module cell_list_mod
    use types_mod
    use config_mod
    use error_mod
    use omp_lib
    implicit none
    private

    !===============================================================================
    ! Module variables
    !===============================================================================
    !> Array of cells in the 3D grid
    type(cell_type), allocatable :: cells(:,:,:)
    
    !> Number of cells in x, y, z dimensions
    integer :: n_cells_x, n_cells_y, n_cells_z
    
    !> Size of each cell
    real(dp) :: cell_size
    
    !> Maximum cutoff distance among all atom pairs
    real(dp) :: max_cutoff
    
    !> Frame number of the last cell update
    integer :: last_cell_update = -1
    
    !> Previous atom coordinates for displacement checking
    real(dp), allocatable :: prev_coords(:,:)
    
    !> Mask indicating which cells need updating (for selective rebuild)
    logical, allocatable :: rebuild_mask(:,:,:)
    
    !> Counter for complete cell list rebuilds
    integer :: num_cell_resets = 0
    
    !> Counter for cell list maintenance updates
    integer :: num_cell_updates = 0
    
    !> Counter for selective cell updates
    integer :: num_selective_updates = 0
    
    !===============================================================================
    ! Cell statistics variables
    !===============================================================================
    !> Number of empty cells in the grid
    integer :: total_empty_cells = 0
    
    !> Maximum number of atoms in any cell
    integer :: max_atoms_in_cell = 0
    
    !> Minimum number of atoms in any non-empty cell
    integer :: min_atoms_in_cell = 0
    
    !> Average number of atoms per non-empty cell
    real(dp) :: avg_atoms_per_cell = 0.0_dp
    
    !> Number of times statistics have been collected
    integer :: stats_collection_count = 0

    !===============================================================================
    ! Public procedures and variables
    !===============================================================================
    public :: initialize_cell_list
    public :: update_cell_list
    public :: get_neighboring_cells
    public :: cleanup_cell_list
    public :: get_cell_grid_dims
    public :: get_cell
    public :: get_cell_statistics
    public :: find_neighbors_within_cutoff
    public :: get_max_displacement
    
    !> Counter for complete cell list rebuilds
    public :: num_cell_resets
    
    !> Counter for cell list maintenance updates
    public :: num_cell_updates
    
    !> Counter for selective cell updates
    public :: num_selective_updates
    
    !> Number of empty cells in the grid
    public :: total_empty_cells
    
    !> Maximum number of atoms in any cell
    public :: max_atoms_in_cell
    
    !> Minimum number of atoms in any non-empty cell
    public :: min_atoms_in_cell
    
    !> Average number of atoms per non-empty cell
    public :: avg_atoms_per_cell
    
    !> Cell size for external reference
    public :: cell_size
    
    !> Grid dimensions for external reference
    public :: n_cells_x, n_cells_y, n_cells_z

contains

    !> @brief Initialize the cell list grid
    !>
    !> Sets up the cell list data structure by dividing the simulation box
    !> into cells of appropriate size based on the maximum cutoff distance.
    !>
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] cutoffs Array of cutoff distances for atom pairs
    !> @param[in] n_atoms Total number of atoms in the system
    subroutine initialize_cell_list(box_length, cutoffs, n_atoms)
        real(dp), intent(in) :: box_length(3)
        real(dp), intent(in) :: cutoffs(:)
        integer, intent(in) :: n_atoms
        integer :: ix, iy, iz, alloc_stat
        
        ! Find maximum cutoff distance
        max_cutoff = maxval(cutoffs)
        
        ! Calculate cell size and grid dimensions
        ! We multiply max_cutoff by CELL_SIZE_FACTOR to get the cell size
        ! Typically CELL_SIZE_FACTOR=1.0, so cells are just large enough to contain
        ! atoms within the cutoff distance
        cell_size = max_cutoff * CELL_SIZE_FACTOR
        
        ! Calculate number of cells in each dimension, ensuring at least 1 cell
        n_cells_x = max(1, floor(box_length(1)/cell_size))
        n_cells_y = max(1, floor(box_length(2)/cell_size))
        n_cells_z = max(1, floor(box_length(3)/cell_size))
        
        ! Adjust cell size to fit box exactly
        cell_size = min(box_length(1)/n_cells_x, &
                       box_length(2)/n_cells_y, &
                       box_length(3)/n_cells_z)

        ! Validate the cell size is at least equal to the maximum cutoff
        if (cell_size < max_cutoff) then
            write(*,*) "WARNING: Cell size smaller than cutoff distance."
            write(*,*) "  Cell size:", cell_size, "Maximum cutoff:", max_cutoff
            write(*,*) "  This may cause missing atom pairs in coordination calculation."
            write(*,*) "  Adjusting cell size to match cutoff."
            
            ! Force cell size to be at least equal to cutoff
            cell_size = max_cutoff
            
            ! Recalculate grid dimensions
            n_cells_x = max(1, floor(box_length(1)/cell_size))
            n_cells_y = max(1, floor(box_length(2)/cell_size))
            n_cells_z = max(1, floor(box_length(3)/cell_size))
        end if

        ! Allocate cell grid
        if (allocated(cells)) deallocate(cells)
        allocate(cells(n_cells_x, n_cells_y, n_cells_z), stat=alloc_stat)
        call check_allocation(alloc_stat, "cell grid")

        ! Allocate rebuild mask for selective rebuilding
        if (allocated(rebuild_mask)) deallocate(rebuild_mask)
        allocate(rebuild_mask(n_cells_x, n_cells_y, n_cells_z), stat=alloc_stat)
        call check_allocation(alloc_stat, "rebuild mask")
        rebuild_mask = .false.

        ! Initialize cells in parallel
        ! We use COLLAPSE(3) to flatten the 3D loop for better load balancing
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    call cells(ix,iy,iz)%init(MAX_ATOMS_PER_CELL)
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        ! Allocate previous coordinates array for displacement tracking
        if (allocated(prev_coords)) deallocate(prev_coords)
        allocate(prev_coords(n_atoms,3), stat=alloc_stat)
        call check_allocation(alloc_stat, "previous coordinates")
        
        ! Reset statistics counters
        num_cell_resets = 0
        num_cell_updates = 0
        num_selective_updates = 0
        total_empty_cells = 0
        max_atoms_in_cell = 0
        min_atoms_in_cell = 0
        avg_atoms_per_cell = 0.0_dp
        stats_collection_count = 0
    end subroutine initialize_cell_list

    !> @brief Update the cell list based on current atom positions
    !>
    !> Updates the cell assignments for atoms based on their current positions.
    !> Supports three update modes:
    !> 1. Complete rebuild (force_update=.true. or force_rebuild=.true.)
    !> 2. Selective rebuild (for atoms that moved significantly)
    !> 3. Skip update (if no atom moved enough)
    !>
    !> @param[in] coords Current atom coordinates
    !> @param[in] atom_types Atom type indices
    !> @param[in] n_atoms Total number of atoms
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] frame Current frame number
    !> @param[in] include_mask Mask indicating which atom types to include
    !> @param[in] force_rebuild Optional flag to force a complete rebuild (default: .false.)
    subroutine update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask, force_rebuild)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        logical, intent(in), optional :: force_rebuild
        
        integer :: ix, iy, iz, i, cell_x, cell_y, cell_z
        real(dp) :: max_displacement, dx, dy, dz, disp_sq
        logical :: force_update, any_selective_updates, do_force_rebuild
        
        ! Initialize variables
        any_selective_updates = .false.
        
        ! Process optional force_rebuild parameter
        do_force_rebuild = .false.
        if (present(force_rebuild)) do_force_rebuild = force_rebuild
        
        ! Check if full update is needed
        force_update = (last_cell_update == -1) .or. do_force_rebuild
        
        ! Custom cell update frequency from setup file
        if (.not. force_update .and. cell_update_freq > 0 .and. &
            mod(frame, cell_update_freq) == 0) then
            force_update = .true.
        end if
        
        if (.not. force_update) then
            ! Check maximum atomic displacement
            max_displacement = get_max_displacement(coords, n_atoms, box_length, include_mask)
            rebuild_mask = .false.
            
            ! If using selective rebuild and displacement is significant but not enough for full rebuild
            if (selective_rebuild .and. max_displacement > REBUILD_THRESHOLD * cell_size * 0.25_dp .and. &
                max_displacement <= REBUILD_THRESHOLD * cell_size) then
                
                ! Mark cells that need updating based on atom movements
                !$OMP PARALLEL DO PRIVATE(dx,dy,dz,disp_sq,cell_x,cell_y,cell_z) SCHEDULE(guided)
                do i = 1, n_atoms
                    ! Skip atoms of excluded types
                    if (.not. include_mask(atom_types(i))) cycle
                    
                    ! Calculate displacement vector
                    dx = coords(i,1) - prev_coords(i,1)
                    dy = coords(i,2) - prev_coords(i,2)
                    dz = coords(i,3) - prev_coords(i,3)
                    
                    ! Apply minimum image convention for periodic boundaries
                    dx = dx - nint(dx/box_length(1)) * box_length(1)
                    dy = dy - nint(dy/box_length(2)) * box_length(2)
                    dz = dz - nint(dz/box_length(3)) * box_length(3)
                    
                    ! Calculate squared displacement
                    disp_sq = dx*dx + dy*dy + dz*dz
                    
                    ! For selective rebuilding, mark cells that need updating
                    ! if atom moved more than half a cell size times the rebuild threshold
                    if (sqrt(disp_sq) > REBUILD_THRESHOLD * cell_size * 0.25_dp) then
                        ! Get current cell indices
                        cell_x = floor(coords(i,1)/cell_size) + 1
                        cell_y = floor(coords(i,2)/cell_size) + 1
                        cell_z = floor(coords(i,3)/cell_size) + 1
                        
                        ! Apply periodic boundary conditions to cell indices
                        cell_x = modulo(cell_x-1, n_cells_x) + 1
                        cell_y = modulo(cell_y-1, n_cells_y) + 1
                        cell_z = modulo(cell_z-1, n_cells_z) + 1
                        
                        ! Mark cell for rebuild - need critical section to avoid race conditions
                        !$OMP CRITICAL(rebuild_mask_update)
                        rebuild_mask(cell_x, cell_y, cell_z) = .true.
                        any_selective_updates = .true.
                        !$OMP END CRITICAL(rebuild_mask_update)
                    end if
                end do
                !$OMP END PARALLEL DO
            end if
            
            ! Force full update if displacement exceeds threshold
            force_update = (max_displacement > REBUILD_THRESHOLD * cell_size)
            
            ! If using selective rebuild but no atoms moved enough, skip update
            if (.not. force_update .and. selective_rebuild .and. .not. any_selective_updates) then
                num_cell_updates = num_cell_updates + 1
                return
            end if
        end if
        
        if (force_update) then
            ! Complete cell reset for full update
            ! Reset all cell counters to zero in parallel
            !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
            do iz = 1, n_cells_z
                do iy = 1, n_cells_y
                    do ix = 1, n_cells_x
                        cells(ix,iy,iz)%n_atoms = 0
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            num_cell_resets = num_cell_resets + 1
            
            ! Add all atoms to cells in parallel
            !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(guided)
            do i = 1, n_atoms
                if (.not. include_mask(atom_types(i))) cycle
                call add_atom_to_cell(i, coords(i,:), box_length)
            end do
            !$OMP END PARALLEL DO
        else if (selective_rebuild) then
            ! Selective cell rebuild
            ! First, mark neighboring cells for update too
            call mark_neighboring_cells()
            
            ! Reset only marked cells in parallel
            !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
            do iz = 1, n_cells_z
                do iy = 1, n_cells_y
                    do ix = 1, n_cells_x
                        if (rebuild_mask(ix,iy,iz)) then
                            cells(ix,iy,iz)%n_atoms = 0
                        end if
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            
            ! Add atoms to marked cells only in parallel
            !$OMP PARALLEL DO PRIVATE(i,cell_x,cell_y,cell_z) SCHEDULE(guided)
            do i = 1, n_atoms
                if (.not. include_mask(atom_types(i))) cycle
                
                ! Get current cell indices
                cell_x = floor(coords(i,1)/cell_size) + 1
                cell_y = floor(coords(i,2)/cell_size) + 1
                cell_z = floor(coords(i,3)/cell_size) + 1
                
                ! Apply periodic boundary conditions
                cell_x = modulo(cell_x-1, n_cells_x) + 1
                cell_y = modulo(cell_y-1, n_cells_y) + 1
                cell_z = modulo(cell_z-1, n_cells_z) + 1
                
                ! Only add atom if its cell needs updating
                if (rebuild_mask(cell_x, cell_y, cell_z)) then
                    call add_atom_to_cell(i, coords(i,:), box_length)
                end if
            end do
            !$OMP END PARALLEL DO
            
            num_selective_updates = num_selective_updates + 1
        end if
        
        ! Update previous coordinates for next displacement check
        prev_coords = coords
        last_cell_update = frame
        
        ! Collect statistics on cell usage
        call collect_cell_statistics()
    end subroutine update_cell_list
    
    !> @brief Calculate maximum atom displacement since last update
    !>
    !> Computes the maximum displacement of any atom relative to its
    !> previous position, applying periodic boundary conditions.
    !>
    !> @param[in] coords Current atom coordinates
    !> @param[in] n_atoms Number of atoms in the system
    !> @param[in] box_length Box dimensions in x, y, z
    !> @param[in] include_mask Mask indicating which atom types to include
    !> @return Maximum displacement of any atom
    function get_max_displacement(coords, n_atoms, box_length, include_mask) result(max_displacement)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        logical, intent(in) :: include_mask(:)
        real(dp) :: max_displacement
        
        integer :: i
        real(dp) :: dx, dy, dz, disp_sq
        
        ! Initialize maximum displacement
        max_displacement = 0.0_dp
        
        ! Calculate maximum atom displacement since last update
        !$OMP PARALLEL DO PRIVATE(dx,dy,dz,disp_sq) REDUCTION(max:max_displacement) SCHEDULE(guided)
        do i = 1, n_atoms
            ! Calculate displacement vector
            dx = coords(i,1) - prev_coords(i,1)
            dy = coords(i,2) - prev_coords(i,2)
            dz = coords(i,3) - prev_coords(i,3)
            
            ! Apply minimum image convention for periodic boundaries
            dx = dx - nint(dx/box_length(1)) * box_length(1)
            dy = dy - nint(dy/box_length(2)) * box_length(2)
            dz = dz - nint(dz/box_length(3)) * box_length(3)
            
            ! Calculate squared displacement
            disp_sq = dx*dx + dy*dy + dz*dz
            
            ! Update maximum displacement
            max_displacement = max(max_displacement, sqrt(disp_sq))
        end do
        !$OMP END PARALLEL DO
        
        return
    end function get_max_displacement
    
    !> @brief Mark neighboring cells for selective rebuild
    !>
    !> When an atom moves significantly, we need to update not only its
    !> current cell but also neighboring cells that might be affected.
    !> This function marks all 26 neighbors of each cell that needs rebuilding.
    subroutine mark_neighboring_cells()
        integer :: ix, iy, iz, dx, dy, dz, nx, ny, nz
        logical :: rebuild_mask_copy(n_cells_x, n_cells_y, n_cells_z)
        
        ! Make a copy to avoid marking new cells during iteration
        rebuild_mask_copy = rebuild_mask
        
        ! Loop through all cells
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    ! If this cell needs rebuilding, mark its neighbors
                    if (rebuild_mask_copy(ix,iy,iz)) then
                        ! Mark all 26 neighboring cells
                        do dz = -1, 1
                            do dy = -1, 1
                                do dx = -1, 1
                                    ! Skip current cell (already marked)
                                    if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                                    
                                    ! Apply periodic boundary conditions
                                    nx = modulo(ix+dx-1, n_cells_x) + 1
                                    ny = modulo(iy+dy-1, n_cells_y) + 1
                                    nz = modulo(iz+dz-1, n_cells_z) + 1
                                    
                                    ! Mark neighboring cell
                                    rebuild_mask(nx, ny, nz) = .true.
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
    end subroutine mark_neighboring_cells

    !> @brief Add an atom to its corresponding cell
    !>
    !> Determines which cell an atom belongs to based on its position
    !> and adds it to that cell. Handles thread-safe updates and dynamic
    !> resizing of cell arrays when needed.
    !>
    !> @param[in] atom_idx Index of the atom to add
    !> @param[in] pos Position of the atom (x,y,z)
    !> @param[in] box_length Box dimensions in x, y, z
    subroutine add_atom_to_cell(atom_idx, pos, box_length)
        integer, intent(in) :: atom_idx
        real(dp), intent(in) :: pos(3)
        real(dp), intent(in) :: box_length(3)
        integer :: ix, iy, iz, current_size, local_n_atoms
        logical :: need_resize
        
        ! Calculate cell indices from atom position
        ix = floor(pos(1)/cell_size) + 1
        iy = floor(pos(2)/cell_size) + 1
        iz = floor(pos(3)/cell_size) + 1
        
        ! Apply periodic boundary conditions
        ix = modulo(ix-1, n_cells_x) + 1
        iy = modulo(iy-1, n_cells_y) + 1
        iz = modulo(iz-1, n_cells_z) + 1
        
        ! Pre-check if resize needed to avoid critical section if possible
        need_resize = .false.
        if (cells(ix,iy,iz)%n_atoms >= size(cells(ix,iy,iz)%atoms)) then
            need_resize = .true.
        end if
        
        ! Only use critical section for resize operation
        if (need_resize) then
            !$OMP CRITICAL(cell_resize)
            current_size = size(cells(ix,iy,iz)%atoms)
            if (cells(ix,iy,iz)%n_atoms >= current_size) then
                ! Double the cell capacity when resizing
                call cells(ix,iy,iz)%resize(2 * current_size)
            end if
            !$OMP END CRITICAL(cell_resize)
        end if
        
        ! Use atomic update for the counter to avoid race conditions
        !$OMP ATOMIC CAPTURE
        cells(ix,iy,iz)%n_atoms = cells(ix,iy,iz)%n_atoms + 1
        local_n_atoms = cells(ix,iy,iz)%n_atoms
        !$OMP END ATOMIC
        
        ! Direct array access outside critical section
        ! Thread-safe because each thread writes to a unique position
        cells(ix,iy,iz)%atoms(local_n_atoms) = atom_idx
    end subroutine add_atom_to_cell

    !> @brief Get list of neighboring cells for a given cell
    !>
    !> Returns the indices of cells that are neighbors of the specified cell,
    !> taking into account periodic boundary conditions.
    !>
    !> @param[in] ix X index of the cell
    !> @param[in] iy Y index of the cell
    !> @param[in] iz Z index of the cell
    !> @param[out] neighbor_cells Array of neighbor cell indices (x,y,z)
    !> @param[out] n_neighbors Number of neighboring cells found
    subroutine get_neighboring_cells(ix, iy, iz, neighbor_cells, n_neighbors)
        integer, intent(in) :: ix, iy, iz
        integer, intent(out) :: neighbor_cells(27,3)
        integer, intent(out) :: n_neighbors
        integer :: dx, dy, dz, nx, ny, nz
        
        n_neighbors = 0
        ! Check all 26 neighboring cells (full 3×3×3 neighborhood excluding center)
        do dz = -1, 1
            do dy = -1, 1
                do dx = -1, 1
                    ! Skip the cell itself (center) as it's handled separately
                    if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                    
                    ! Apply periodic boundary conditions
                    nx = modulo(ix+dx-1, n_cells_x) + 1
                    ny = modulo(iy+dy-1, n_cells_y) + 1
                    nz = modulo(iz+dz-1, n_cells_z) + 1
                    
                    n_neighbors = n_neighbors + 1
                    neighbor_cells(n_neighbors,1) = nx
                    neighbor_cells(n_neighbors,2) = ny
                    neighbor_cells(n_neighbors,3) = nz
                end do
            end do
        end do
    end subroutine get_neighboring_cells

    !> @brief Get all cells within a certain radius of a position
    !>
    !> Returns the indices of all cells that could contain atoms within
    !> the specified cutoff distance from the given position.
    !>
    !> @param[in] pos Position to find cells around
    !> @param[in] radius Distance radius to consider
    !> @param[in] box_length Box dimensions
    !> @param[out] neighbor_cells Array of cell indices (x,y,z)
    !> @param[out] n_cells Number of cells found
    subroutine get_cells_within_radius(pos, radius, box_length, neighbor_cells, n_cells)
        real(dp), intent(in) :: pos(3)
        real(dp), intent(in) :: radius
        real(dp), intent(in) :: box_length(3)
        integer, intent(out) :: neighbor_cells(:,:)
        integer, intent(out) :: n_cells
        
        integer :: center_x, center_y, center_z
        integer :: min_x, max_x, min_y, max_y, min_z, max_z
        integer :: nx, ny, nz, ix, iy, iz
        integer :: cells_per_radius
        
        ! Calculate the central cell containing the position
        center_x = floor(pos(1)/cell_size) + 1
        center_y = floor(pos(2)/cell_size) + 1
        center_z = floor(pos(3)/cell_size) + 1
        
        ! Apply periodic boundary conditions
        center_x = modulo(center_x-1, n_cells_x) + 1
        center_y = modulo(center_y-1, n_cells_y) + 1
        center_z = modulo(center_z-1, n_cells_z) + 1
        
        ! Calculate how many cells to go in each direction
        cells_per_radius = ceiling(radius / cell_size)
        
        ! Calculate min/max cell indices with periodic boundaries
        min_x = center_x - cells_per_radius
        max_x = center_x + cells_per_radius
        min_y = center_y - cells_per_radius
        max_y = center_y + cells_per_radius
        min_z = center_z - cells_per_radius
        max_z = center_z + cells_per_radius
        
        ! Initialize count
        n_cells = 0
        
        ! Collect all cells in the range
        do iz = min_z, max_z
            do iy = min_y, max_y
                do ix = min_x, max_x
                    ! Apply periodic boundary conditions
                    nx = modulo(ix-1, n_cells_x) + 1
                    ny = modulo(iy-1, n_cells_y) + 1
                    nz = modulo(iz-1, n_cells_z) + 1
                    
                    ! Add cell to list if we have space
                    if (n_cells < size(neighbor_cells, 1)) then
                        n_cells = n_cells + 1
                        neighbor_cells(n_cells, 1) = nx
                        neighbor_cells(n_cells, 2) = ny
                        neighbor_cells(n_cells, 3) = nz
                    else
                        ! Not enough space, return what we have
                        return
                    end if
                end do
            end do
        end do
    end subroutine get_cells_within_radius
    
!> @brief Find all atoms within a cutoff distance of a given atom
!>
!> Uses the cell list to efficiently identify all atoms within the
!> specified cutoff distance of the target atom, applying periodic
!> boundary conditions. Primarily used for building Verlet lists.
!>
!> @param[in] atom_idx Index of the atom to find neighbors for
!> @param[in] coords Atom coordinates
!> @param[in] atom_types Atom type indices
!> @param[in] box_length Box dimensions
!> @param[in] cutoff_distance Cutoff distance for neighbor finding
!> @param[out] neighbors Array to store neighbor atom indices
!> @param[out] pair_types Array to store pair type indices for each neighbor
!> @param[out] n_neighbors Number of neighbors found
!> @param[in] max_neighbors Maximum number of neighbors to return
!> @param[in] pairs Type pair definitions
!> @param[in] n_pairs Number of pair types
subroutine find_neighbors_within_cutoff(atom_idx, coords, atom_types, box_length, &
                                      cutoff_distance, neighbors, pair_types, &
                                      n_neighbors, max_neighbors, pairs, n_pairs)
    integer, intent(in) :: atom_idx
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: atom_types(:)
    real(dp), intent(in) :: box_length(3)
    real(dp), intent(in) :: cutoff_distance
    integer, intent(out) :: neighbors(:)
    integer, intent(out) :: pair_types(:)
    integer, intent(out) :: n_neighbors
    integer, intent(in) :: max_neighbors
    type(pair_type), intent(in) :: pairs(:)
    integer, intent(in) :: n_pairs
    
    integer :: neighbor_cells(125,3), n_cells, cell_idx
    integer :: cell_ix, cell_iy, cell_iz
    integer :: atom_j, type_i, type_j, pair_idx, atom_loop_idx
    real(dp) :: pos(3), dx, dy, dz, r_sq, pair_cutoff_sq
    
    ! Initialize
    n_neighbors = 0
    pos = coords(atom_idx,:)
    type_i = atom_types(atom_idx)
    
    ! Get all cells within the cutoff radius (using the Verlet cutoff + skin)
    call get_cells_within_radius(pos, cutoff_distance, box_length, neighbor_cells, n_cells)
    
    ! Loop through all cells
    do cell_idx = 1, n_cells
        ! Get current cell indices
        cell_ix = neighbor_cells(cell_idx, 1)
        cell_iy = neighbor_cells(cell_idx, 2)
        cell_iz = neighbor_cells(cell_idx, 3)
        
        ! Loop through atoms in the cell
        do atom_loop_idx = 1, cells(cell_ix, cell_iy, cell_iz)%n_atoms
            atom_j = cells(cell_ix, cell_iy, cell_iz)%atoms(atom_loop_idx)
            
            ! Skip self
            if (atom_j == atom_idx) cycle
            
            ! Get atom types and find pair type index
            type_j = atom_types(atom_j)
            call get_pair_index(type_i, type_j, pairs, n_pairs, pair_idx)
            
            ! Skip if no valid pair type found
            if (pair_idx <= 0) cycle
            
            ! Calculate distance
            dx = coords(atom_j,1) - pos(1)
            dy = coords(atom_j,2) - pos(2)
            dz = coords(atom_j,3) - pos(3)
            
            ! Apply minimum image convention
            dx = dx - nint(dx/box_length(1)) * box_length(1)
            dy = dy - nint(dy/box_length(2)) * box_length(2)
            dz = dz - nint(dz/box_length(3)) * box_length(3)
            
            ! Calculate squared distance
            r_sq = dx*dx + dy*dy + dz*dz
            
            ! Use the pair-specific cutoff instead of global cutoff
            pair_cutoff_sq = pairs(pair_idx)%cutoff_sq
            
            ! Check if within the pair's specific cutoff plus skin distance
            ! The cutoff_distance parameter includes the skin distance
            if (r_sq <= cutoff_distance * cutoff_distance) then
                ! Check if we have space and add neighbor
                if (n_neighbors < max_neighbors) then
                    n_neighbors = n_neighbors + 1
                    neighbors(n_neighbors) = atom_j
                    pair_types(n_neighbors) = pair_idx
                else
                    ! Buffer is full, exit
                    return
                end if
            end if
        end do
    end do
end subroutine find_neighbors_within_cutoff
    
    !> @brief Get the pair index for two atom types
    !>
    !> Returns the index in the pairs array that corresponds to the
    !> interaction between two atoms with the specified type indices.
    !>
    !> @param[in] type1 Type index of first atom
    !> @param[in] type2 Type index of second atom
    !> @param[in] pairs Array of pair type definitions
    !> @param[in] n_pairs Number of pair types
    !> @param[out] pair_idx Index of the pair type (-1 if not found)
    subroutine get_pair_index(type1, type2, pairs, n_pairs, pair_idx)
        integer, intent(in) :: type1, type2
        type(pair_type), intent(in) :: pairs(:)
        integer, intent(in) :: n_pairs
        integer, intent(out) :: pair_idx
        
        integer :: i, t1, t2
        
        ! Ensure types are ordered consistently (smaller first)
        t1 = min(type1, type2)
        t2 = max(type1, type2)
        
        ! Search for matching pair
        do i = 1, n_pairs
            if (pairs(i)%type1_id == t1 .and. pairs(i)%type2_id == t2) then
                pair_idx = i
                return
            end if
        end do
        
        ! No matching pair found
        pair_idx = -1
    end subroutine get_pair_index

    !> @brief Get the dimensions of the cell grid
    !>
    !> Returns the number of cells in each dimension.
    !>
    !> @param[out] nx Number of cells in x direction
    !> @param[out] ny Number of cells in y direction
    !> @param[out] nz Number of cells in z direction
    subroutine get_cell_grid_dims(nx, ny, nz)
        integer, intent(out) :: nx, ny, nz
        nx = n_cells_x
        ny = n_cells_y
        nz = n_cells_z
    end subroutine get_cell_grid_dims

    !> @brief Get a specific cell from the grid
    !>
    !> Returns a copy of the cell at the specified indices.
    !>
    !> @param[in] ix X index of the cell
    !> @param[in] iy Y index of the cell
    !> @param[in] iz Z index of the cell
    !> @param[out] cell The requested cell
    subroutine get_cell(ix, iy, iz, cell)
        integer, intent(in) :: ix, iy, iz
        type(cell_type), intent(out) :: cell
        cell = cells(ix, iy, iz)
    end subroutine get_cell

    !> @brief Calculate unique cell index from 3D indices
    !>
    !> Helper function to create a unique 1D index for each cell
    !> from its 3D grid coordinates. Used for cell pair ordering.
    !>
    !> @param[in] ix X index of the cell
    !> @param[in] iy Y index of the cell
    !> @param[in] iz Z index of the cell
    !> @return Unique 1D index for the cell
    function cell_index(ix, iy, iz) result(idx)
        integer, intent(in) :: ix, iy, iz
        integer :: idx
        
        idx = ix + (iy-1)*n_cells_x + (iz-1)*n_cells_x*n_cells_y
    end function cell_index
    
    !> @brief Collect statistics on cell usage
    !>
    !> Gathers information about cell occupancy, including the number of
    !> empty cells, and the minimum, maximum, and average atoms per cell.
    subroutine collect_cell_statistics()
        integer :: ix, iy, iz, empty_count, total_cells, total_atoms
        integer :: local_min_atoms, local_max_atoms
        
        ! Initialize counters
        empty_count = 0
        total_cells = n_cells_x * n_cells_y * n_cells_z
        total_atoms = 0
        local_min_atoms = huge(local_min_atoms)
        local_max_atoms = 0
        
        ! Count statistics
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    ! Track empty cells
                    if (cells(ix,iy,iz)%n_atoms == 0) then
                        empty_count = empty_count + 1
                    end if
                    
                    ! Track min/max atoms in a cell
                    local_max_atoms = max(local_max_atoms, cells(ix,iy,iz)%n_atoms)
                    if (cells(ix,iy,iz)%n_atoms > 0) then
                        local_min_atoms = min(local_min_atoms, cells(ix,iy,iz)%n_atoms)
                    end if
                    
                    ! Track total atoms for average calculation
                    total_atoms = total_atoms + cells(ix,iy,iz)%n_atoms
                end do
            end do
        end do
        
        ! Update global statistics (running average)
        stats_collection_count = stats_collection_count + 1
        
        ! Update max and min (global maxima and minima)
        if (stats_collection_count == 1) then
            ! First collection, initialize values
            min_atoms_in_cell = local_min_atoms
            max_atoms_in_cell = local_max_atoms
            total_empty_cells = empty_count
            avg_atoms_per_cell = real(total_atoms, dp) / max(1, total_cells - empty_count)
        else
            ! Update running statistics
            min_atoms_in_cell = min(min_atoms_in_cell, local_min_atoms)
            max_atoms_in_cell = max(max_atoms_in_cell, local_max_atoms)
            
            ! Update running average of empty cells and atoms per cell
            total_empty_cells = nint((total_empty_cells * (stats_collection_count - 1) + &
                                    empty_count) / real(stats_collection_count, dp))
            
            ! Calculate average atoms per non-empty cell
            avg_atoms_per_cell = ((avg_atoms_per_cell * (stats_collection_count - 1)) + &
                                real(total_atoms, dp) / max(1, total_cells - empty_count)) / &
                                real(stats_collection_count, dp)
        end if
    end subroutine collect_cell_statistics
    
    !> @brief Get statistics on cell usage
    !>
    !> Returns an array of statistics about cell usage.
    !>
    !> @return Array of statistics:
    !>     (1) Number of empty cells
    !>     (2) Average atoms per non-empty cell
    !>     (3) Minimum atoms in any non-empty cell
    !>     (4) Maximum atoms in any cell
    function get_cell_statistics() result(stats)
        real(dp) :: stats(4)
        stats(1) = real(total_empty_cells, dp)
        stats(2) = avg_atoms_per_cell
        stats(3) = real(min_atoms_in_cell, dp)
        stats(4) = real(max_atoms_in_cell, dp)
    end function get_cell_statistics

    !> @brief Clean up cell list resources
    !>
    !> Deallocates all memory used by the cell list, including the cell grid,
    !> previous coordinates array, and rebuild mask.
    subroutine cleanup_cell_list()
        integer :: ix, iy, iz
        
        if (allocated(cells)) then
            ! Clean up each cell in parallel
            !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
            do iz = 1, n_cells_z
                do iy = 1, n_cells_y
                    do ix = 1, n_cells_x
                        call cells(ix,iy,iz)%cleanup()
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            deallocate(cells)
        end if
        
        if (allocated(prev_coords)) deallocate(prev_coords)
        if (allocated(rebuild_mask)) deallocate(rebuild_mask)
    end subroutine cleanup_cell_list

end module cell_list_mod
