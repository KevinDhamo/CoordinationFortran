!> @brief Module defining fundamental data types for coordination analysis
!>
!> This module provides the core data type definitions used throughout
!> the coordination analysis program, including atom type information,
!> pair interaction specifications, and cell data structures for spatial
!> partitioning algorithms. It also exports the high-precision real kind (dp).
module types_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    
    ! Public type definitions
    public :: atom_type_info, pair_type, cell_type, verlet_list_type
    public :: dp

    !===============================================================================
    ! Type definitions
    !===============================================================================
    
    !> @brief Information about atom types in the system
    !>
    !> Contains properties of atom types including numeric ID, element name,
    !> mass, and a flag indicating whether to include this type in analysis.
    type :: atom_type_info
        integer :: type_id           !< Numeric type ID from data file
        character(len=8) :: name     !< Element name from data file comment
        real(dp) :: mass             !< Mass from data file
        logical :: include = .true.  !< Whether to include in analysis
    end type atom_type_info
    
    !> @brief Pair interaction specification
    !>
    !> Defines the interaction between two atom types, including their IDs
    !> and the cutoff distance for considering them as neighbors.
    type :: pair_type
        integer :: type1_id, type2_id  !< Numeric type IDs for the pair
        real(dp) :: cutoff             !< Cutoff distance for coordination
        real(dp) :: cutoff_sq          !< Squared cutoff (for optimization)
    end type pair_type

    !> @brief Cell data structure for spatial partitioning
    !>
    !> Represents a spatial cell in the cell list algorithm used for
    !> efficient neighbor finding. Contains a list of atoms in the cell.
    type :: cell_type
        integer :: n_atoms = 0                !< Number of atoms in cell
        integer, allocatable :: atoms(:)      !< Indices of atoms in this cell
    contains
        !> Initialize the cell with a specified size
        procedure :: init => initialize_cell
        
        !> Clean up cell resources
        procedure :: cleanup => cleanup_cell
        
        !> Resize the atom array when needed
        procedure :: resize => resize_cell_array
    end type cell_type
    
    !> @brief Verlet neighbor list for an atom
    !>
    !> Contains a list of neighboring atoms within a cutoff plus skin distance.
    !> Used for efficient coordination number calculations by avoiding the need
    !> to search for neighbors at every step.
    type :: verlet_list_type
        integer :: n_neighbors = 0               !< Number of neighbors in list
        integer, allocatable :: neighbors(:)     !< Array of neighbor atom indices
        integer, allocatable :: pair_types(:)    !< Pair type index for each neighbor
    contains
        !> Initialize the Verlet list with a specified capacity
        procedure :: init => initialize_verlet_list
        
        !> Clean up Verlet list resources
        procedure :: cleanup => cleanup_verlet_list
        
        !> Resize the neighbors array when needed
        procedure :: resize => resize_verlet_list
    end type verlet_list_type

contains
    
    !> @brief Initialize a cell with specified initial capacity
    !>
    !> Allocates memory for the atom indices array with the specified
    !> initial size and sets the atom count to zero.
    !>
    !> @param[in,out] this Cell to initialize
    !> @param[in] initial_size Initial allocation size for atom array
    subroutine initialize_cell(this, initial_size)
        class(cell_type), intent(inout) :: this
        integer, intent(in) :: initial_size
        integer :: alloc_stat
        
        ! Clean up existing allocation if present
        if (allocated(this%atoms)) deallocate(this%atoms)
        
        ! Allocate atom array with specified initial size
        allocate(this%atoms(initial_size), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Failed to allocate cell array"
        
        ! Initialize atom count
        this%n_atoms = 0
    end subroutine initialize_cell
    
    !> @brief Clean up cell resources
    !>
    !> Deallocates the atom indices array and resets the atom count.
    !>
    !> @param[in,out] this Cell to clean up
    subroutine cleanup_cell(this)
        class(cell_type), intent(inout) :: this
        
        ! Deallocate atom array if allocated
        if (allocated(this%atoms)) deallocate(this%atoms)
        
        ! Reset atom count
        this%n_atoms = 0
    end subroutine cleanup_cell
    
    !> @brief Resize the cell's atom array to accommodate more atoms
    !>
    !> Creates a new, larger array, copies existing data, and replaces
    !> the cell's atom array with the new one.
    !>
    !> @param[in,out] this Cell whose array to resize
    !> @param[in] new_size New size for the atom array
    subroutine resize_cell_array(this, new_size)
        class(cell_type), intent(inout) :: this
        integer, intent(in) :: new_size
        integer, allocatable :: temp(:)
        integer :: alloc_stat
        
        ! Allocate temporary array with new size
        allocate(temp(new_size), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Failed to allocate temporary array"
        
        ! Copy existing data to temporary array
        temp(1:this%n_atoms) = this%atoms(1:this%n_atoms)
        
        ! Replace old array with new one using move_alloc
        ! (avoids copy and automatically deallocates old array)
        call move_alloc(from=temp, to=this%atoms)
    end subroutine resize_cell_array

    !> @brief Initialize a Verlet list with specified initial capacity
    !>
    !> Allocates memory for the neighbor indices array with the specified
    !> initial size and sets the neighbor count to zero.
    !>
    !> @param[in,out] this Verlet list to initialize
    !> @param[in] initial_size Initial allocation size for neighbor arrays
    subroutine initialize_verlet_list(this, initial_size)
        class(verlet_list_type), intent(inout) :: this
        integer, intent(in) :: initial_size
        integer :: alloc_stat
        
        ! Clean up existing allocations if present
        if (allocated(this%neighbors)) deallocate(this%neighbors)
        if (allocated(this%pair_types)) deallocate(this%pair_types)
        
        ! Allocate arrays with specified initial size
        allocate(this%neighbors(initial_size), this%pair_types(initial_size), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Failed to allocate Verlet list arrays"
        
        ! Initialize neighbor count
        this%n_neighbors = 0
    end subroutine initialize_verlet_list
    
    !> @brief Clean up Verlet list resources
    !>
    !> Deallocates the neighbor arrays and resets the neighbor count.
    !>
    !> @param[in,out] this Verlet list to clean up
    subroutine cleanup_verlet_list(this)
        class(verlet_list_type), intent(inout) :: this
        
        ! Deallocate arrays if allocated
        if (allocated(this%neighbors)) deallocate(this%neighbors)
        if (allocated(this%pair_types)) deallocate(this%pair_types)
        
        ! Reset neighbor count
        this%n_neighbors = 0
    end subroutine cleanup_verlet_list
    
    !> @brief Resize the Verlet list's neighbor arrays to accommodate more neighbors
    !>
    !> Creates new, larger arrays, copies existing data, and replaces
    !> the Verlet list's arrays with the new ones.
    !>
    !> @param[in,out] this Verlet list whose arrays to resize
    !> @param[in] new_size New size for the neighbor arrays
    subroutine resize_verlet_list(this, new_size)
        class(verlet_list_type), intent(inout) :: this
        integer, intent(in) :: new_size
        integer, allocatable :: temp_neighbors(:), temp_pair_types(:)
        integer :: alloc_stat
        
        ! Allocate temporary arrays with new size
        allocate(temp_neighbors(new_size), temp_pair_types(new_size), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Failed to allocate temporary Verlet list arrays"
        
        ! Copy existing data to temporary arrays
        temp_neighbors(1:this%n_neighbors) = this%neighbors(1:this%n_neighbors)
        temp_pair_types(1:this%n_neighbors) = this%pair_types(1:this%n_neighbors)
        
        ! Replace old arrays with new ones using move_alloc
        call move_alloc(from=temp_neighbors, to=this%neighbors)
        call move_alloc(from=temp_pair_types, to=this%pair_types)
    end subroutine resize_verlet_list

end module types_mod
