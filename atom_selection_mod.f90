!> @brief Module for selecting atoms based on ID ranges
!>
!> This module provides functionality to select atoms based on a selection string
!> that specifies atom ID ranges (e.g., "1-100,200,300-400").
!> If no selection is active, all atoms are considered selected.
module atom_selection_mod
    use types_mod
    use config_mod
    use error_mod
    implicit none
    private

    ! Public procedures
    public :: initialize_atom_selection
    public :: parse_atom_selection
    public :: is_atom_selected
    public :: get_selected_atom_ids
    public :: get_selected_atom_count
    public :: cleanup_atom_selection
    
    ! Public variables
    public :: selection_active

    ! Module variables
    integer, allocatable :: selected_atom_ids(:)  !< Array of selected atom IDs
    integer :: n_selected_atoms = 0             !< Count of selected atoms
    logical :: selection_active = .false.       !< Flag indicating if selection is active
    
    ! Track if selection message has been shown
    logical :: selection_message_shown = .false.

contains
    !> @brief Initialize the atom selection
    !>
    !> Initializes the atom selection system for n_atoms atoms.
    !>
    !> @param[in] n_atoms Number of atoms in the system
    subroutine initialize_atom_selection(n_atoms)
        integer, intent(in) :: n_atoms
        integer :: alloc_stat
        
        ! Deallocate if already allocated
        if (allocated(selected_atom_ids)) deallocate(selected_atom_ids)
        
        ! Initially allocate enough space for all atoms (will be resized after parsing)
        allocate(selected_atom_ids(n_atoms), stat=alloc_stat)
        call check_allocation(alloc_stat, "selected_atom_ids")
        
        ! Initialize with zeros
        selected_atom_ids = 0
        n_selected_atoms = 0
        selection_active = .false.
        selection_message_shown = .false.  ! Reset message flag
    end subroutine initialize_atom_selection

    !> @brief Parse a selection string and update the atom selection
    !>
    !> Parses a selection string that specifies atom ID ranges, e.g., "1-100,200,300-400"
    !> and builds the list of selected atom IDs. If the selection string is empty
    !> or "all", all atoms are considered selected.
    !>
    !> @param[in] selection_string Selection string specifying atom ID ranges
    !> @param[in] n_atoms Number of atoms in the system
    subroutine parse_atom_selection(selection_string, n_atoms)
        character(len=*), intent(in) :: selection_string
        integer, intent(in) :: n_atoms
        integer :: i, j, comma_pos, dash_pos, start_id, end_id, read_stat
        character(len=256) :: remaining, token
        logical, allocatable :: temp_mask(:)
        integer :: max_selected_atoms
        
        ! If selection string is empty or "all", all atoms are selected
        if (len_trim(selection_string) == 0 .or. &
            trim(selection_string) == 'all' .or. &
            trim(selection_string) == 'ALL') then
            
            selection_active = .false.
            n_selected_atoms = n_atoms
            
            ! Fill with sequential IDs from 1 to n_atoms
            do i = 1, n_atoms
                selected_atom_ids(i) = i
            end do
            
            if (VERBOSE .and. .not. selection_message_shown) then
                write(*,'(A)') " All atoms selected"
                selection_message_shown = .true.
            end if
            return
        end if
        
        ! Use a temporary mask to efficiently track selections
        allocate(temp_mask(n_atoms))
        temp_mask = .false.
        
        ! Process selection string - format: "1-100,200,300-400"
        remaining = adjustl(selection_string)
        selection_active = .true.
        
        do while (len_trim(remaining) > 0)
            ! Find next comma
            comma_pos = index(remaining, ',')
            
            if (comma_pos > 0) then
                token = adjustl(remaining(:comma_pos-1))
                remaining = adjustl(remaining(comma_pos+1:))
            else
                token = adjustl(remaining)
                remaining = ''
            end if
            
            ! Check if token contains a range (dash)
            dash_pos = index(token, '-')
            
            if (dash_pos > 0) then
                ! Process range (e.g., 1-100)
                read(token(:dash_pos-1), *, iostat=read_stat) start_id
                if (read_stat /= 0) cycle
                
                read(token(dash_pos+1:), *, iostat=read_stat) end_id
                if (read_stat /= 0) cycle
                
                ! Validate range
                start_id = max(1, start_id)
                end_id = min(n_atoms, end_id)
                
                ! Mark atoms in range as selected
                do i = start_id, end_id
                    temp_mask(i) = .true.
                end do
            else
                ! Process single atom ID
                read(token, *, iostat=read_stat) i
                if (read_stat /= 0) cycle
                
                ! Validate atom ID
                if (i >= 1 .and. i <= n_atoms) then
                    temp_mask(i) = .true.
                end if
            end if
        end do
        
        ! Count selected atoms and build selected_atom_ids array
        n_selected_atoms = 0
        do i = 1, n_atoms
            if (temp_mask(i)) then
                n_selected_atoms = n_selected_atoms + 1
                selected_atom_ids(n_selected_atoms) = i
            end if
        end do
        
        if (VERBOSE .and. .not. selection_message_shown) then
            write(*,'(A,I0,A,I0,A)') " Selected ", n_selected_atoms, " atoms out of ", n_atoms, " total"
            selection_message_shown = .true.
        end if
        
        ! If no valid atoms were selected, select all as a fallback
        if (n_selected_atoms == 0) then
            if (VERBOSE .and. .not. selection_message_shown) then
                write(*,'(A)') " Warning: No valid atoms selected, defaulting to all atoms"
            end if
            
            ! Reset to include all atoms
            n_selected_atoms = n_atoms
            do i = 1, n_atoms
                selected_atom_ids(i) = i
            end do
            selection_active = .false.
        else
            ! Resize the selected_atom_ids array to save memory
            max_selected_atoms = n_selected_atoms
            call resize_selected_atoms(max_selected_atoms)
        end if
        
        ! Clean up
        deallocate(temp_mask)
    end subroutine parse_atom_selection

    !> @brief Resize the selected atoms array
    !>
    !> Helper function to resize the selected_atom_ids array to the specified size.
    !>
    !> @param[in] new_size New size for the array
    subroutine resize_selected_atoms(new_size)
        integer, intent(in) :: new_size
        integer, allocatable :: temp(:)
        integer :: alloc_stat
        
        ! Create temporary array
        allocate(temp(new_size), stat=alloc_stat)
        call check_allocation(alloc_stat, "temp array for selected_atom_ids resize")
        
        ! Copy data to temporary array (only up to the minimum of old and new sizes)
        temp(1:min(n_selected_atoms, new_size)) = selected_atom_ids(1:min(n_selected_atoms, new_size))
        
        ! Replace original array with resized version
        call move_alloc(from=temp, to=selected_atom_ids)
    end subroutine resize_selected_atoms

    !> @brief Check if an atom is selected
    !>
    !> @param[in] atom_id Atom ID to check
    !> @return .true. if the atom is selected, .false. otherwise
    function is_atom_selected(atom_id) result(selected)
        integer, intent(in) :: atom_id
        logical :: selected
        integer :: i
        
        ! If no selection is active, all atoms are selected
        if (.not. selection_active) then
            selected = .true.
            return
        end if
        
        ! Check if atom_id is in the selected_atom_ids array
        ! This is a linear search, which is not ideal for large selections
        selected = .false.
        do i = 1, n_selected_atoms
            if (selected_atom_ids(i) == atom_id) then
                selected = .true.
                return
            end if
        end do
    end function is_atom_selected

    !> @brief Get the list of selected atom IDs
    !>
    !> Returns a copy of the internal array of selected atom IDs.
    !>
    !> @param[out] atom_ids Array to store the selected atom IDs
    !> @param[out] count Number of selected atoms
    subroutine get_selected_atom_ids(atom_ids, count)
        integer, intent(out) :: atom_ids(:)
        integer, intent(out) :: count
        integer :: copy_size
        
        ! Determine how many elements to copy
        copy_size = min(size(atom_ids), n_selected_atoms)
        
        ! Copy the data
        if (copy_size > 0) then
            atom_ids(1:copy_size) = selected_atom_ids(1:copy_size)
        end if
        
        count = n_selected_atoms
    end subroutine get_selected_atom_ids

    !> @brief Get the count of selected atoms
    !>
    !> @return Number of selected atoms
    function get_selected_atom_count() result(count)
        integer :: count
        count = n_selected_atoms
    end function get_selected_atom_count

    !> @brief Clean up atom selection resources
    subroutine cleanup_atom_selection()
        if (allocated(selected_atom_ids)) deallocate(selected_atom_ids)
        n_selected_atoms = 0
        selection_active = .false.
        selection_message_shown = .false.
    end subroutine cleanup_atom_selection

end module atom_selection_mod