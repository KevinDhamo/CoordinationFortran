module ram_cache_mod
    use types_mod
    use config_mod, only: VERBOSE
    implicit none
    private

    public :: coordination_cache, create_cache, destroy_cache

    type :: cache_entry
        integer :: frame_number = -1
        integer, allocatable :: coordination_numbers(:,:)
        real(dp), allocatable :: coords(:,:)
        integer, allocatable :: atom_types(:)
        character(len=8), allocatable :: elements(:)
    end type cache_entry

    type :: coordination_cache
        type(cache_entry), allocatable :: entries(:)
        integer :: max_entries = 10
        integer :: current_entries = 0
        integer :: replacement_strategy = 1 ! LRU
    contains
        procedure :: add_frame
        procedure :: get_frame
        procedure :: clear
    end type coordination_cache

contains
    function create_cache(max_entries, n_atoms, n_pairs) result(cache)
        integer, intent(in) :: max_entries
        integer, intent(in) :: n_atoms
        integer, intent(in) :: n_pairs
        type(coordination_cache) :: cache
        integer :: i, alloc_stat

        allocate(cache%entries(max_entries), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate cache entries", ERR_ALLOCATION)
            return
        end if

        ! Pre-allocate arrays for each entry
        do i = 1, max_entries
            allocate(cache%entries(i)%coordination_numbers(n_atoms, n_pairs), stat=alloc_stat)
            allocate(cache%entries(i)%coords(n_atoms, 3), stat=alloc_stat)
            allocate(cache%entries(i)%atom_types(n_atoms), stat=alloc_stat)
            allocate(cache%entries(i)%elements(n_atoms), stat=alloc_stat)
        end do

        cache%max_entries = max_entries
        cache%current_entries = 0
    end function create_cache

    subroutine add_frame(this, frame_number, coord_numbers, coords, atom_types, elements)
        class(coordination_cache), intent(inout) :: this
        integer, intent(in) :: frame_number
        integer, intent(in) :: coord_numbers(:,:)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        character(len=*), intent(in) :: elements(:)
        integer :: replace_index, i

        ! Check if frame already exists
        do i = 1, this%current_entries
            if (this%entries(i)%frame_number == frame_number) then
                ! Update existing entry
                this%entries(i)%coordination_numbers = coord_numbers
                this%entries(i)%coords = coords
                this%entries(i)%atom_types = atom_types
                this%entries(i)%elements = elements
                return
            end if
        end do

        ! Determine replacement index
        if (this%current_entries < this%max_entries) then
            this%current_entries = this%current_entries + 1
            replace_index = this%current_entries
        else
            ! Simple round-robin replacement
            replace_index = mod(this%current_entries, this%max_entries) + 1
        end if

        ! Replace entry
        this%entries(replace_index)%frame_number = frame_number
        this%entries(replace_index)%coordination_numbers = coord_numbers
        this%entries(replace_index)%coords = coords
        this%entries(replace_index)%atom_types = atom_types
        this%entries(replace_index)%elements = elements
    end subroutine add_frame

    function get_frame(this, frame_number) result(frame_data)
        class(coordination_cache), intent(in) :: this
        integer, intent(in) :: frame_number
        type(cache_entry) :: frame_data
        integer :: i

        frame_data%frame_number = -1
        do i = 1, this%current_entries
            if (this%entries(i)%frame_number == frame_number) then
                frame_data = this%entries(i)
                return
            end if
        end do
    end function get_frame

    subroutine clear(this)
        class(coordination_cache), intent(inout) :: this
        this%current_entries = 0
    end subroutine clear
end module ram_cache_mod