!> @brief Module for reading molecular dynamics trajectory files
!>
!> This module provides functionality to read molecular dynamics trajectory
!> files in different formats (LAMMPS dump and XYZ). It supports counting
!> frames, seeking to specific frames, and reading frame data with proper
!> handling of periodic boundary conditions.
module trajectory_io_mod
    use types_mod
    use config_mod
    use error_mod
    implicit none
    private

    !===============================================================================
    ! Public procedures
    !===============================================================================
    public :: init_trajectory_reader
    public :: read_trajectory_frame
    public :: count_trajectory_frames
    public :: seek_to_frame
    public :: cleanup_trajectory_reader
    public :: set_user_header_line

    !===============================================================================
    ! Module variables
    !===============================================================================
    !> File unit for input trajectory file
    integer :: input_unit = -1
    
    !> Flag indicating if input file is open
    logical :: input_file_open = .false.
    
    !===============================================================================
    ! Column mapping type definition
    !===============================================================================
    !> Type for mapping trajectory columns to data fields
    type :: column_map_type
        integer :: id = -1      !< Column index for atom ID
        integer :: type = -1    !< Column index for atom type
        integer :: x = -1       !< Column index for x coordinate
        integer :: y = -1       !< Column index for y coordinate
        integer :: z = -1       !< Column index for z coordinate
        integer :: max_cols = -1 !< Maximum column index needed
    end type
    
    !> Current column mapping
    type(column_map_type), save :: col_map
    
    !===============================================================================
    ! Frame location tracking
    !===============================================================================
    !> Array of file offsets for each frame (for fast seeking)
    integer, allocatable :: frame_offsets(:)
    
    !> Flag indicating if frame offsets are available
    logical :: have_frame_offsets = .false.
    
    !> Buffer for reading trajectory values
    real(dp), allocatable :: values_buffer(:)
    
    !===============================================================================
    ! File format types
    !===============================================================================
    !> LAMMPS dump file format identifier
    integer, parameter :: FORMAT_LAMMPS = 1
    
    !> XYZ file format identifier
    integer, parameter :: FORMAT_XYZ = 2
    
    !> Current file format
    integer :: file_format = FORMAT_LAMMPS

    !===============================================================================
    ! Header format options
    !===============================================================================
    !> User-provided header line for custom column mapping
    character(len=512) :: user_header_line = ""
    
    !> Flag to use user-provided header
    logical :: use_user_header = .false.
    
    !> Flag to track if header info has been shown
    logical :: header_shown = .false.
    
    !> Flag to track if column info has been shown
    logical :: shown_column_info = .false.

contains

    !> @brief Initialize trajectory file reader
    !>
    !> Opens the trajectory file, determines its format based on extension,
    !> and sets up initial state for reading frames.
    !>
    !> @param[in] filename Path to trajectory file
    subroutine init_trajectory_reader(filename)
        character(len=*), intent(in) :: filename
        integer :: io_stat
        character(len=10) :: file_ext
        
        ! Reset tracking variables when initializing reader
        header_shown = .false.
        shown_column_info = .false.
        
        ! Check file extension to determine format
        file_ext = get_file_extension(filename)
        
        if (trim(file_ext) == "xyz") then
            file_format = FORMAT_XYZ
        else
            file_format = FORMAT_LAMMPS
        end if
        
        ! Open trajectory file
        open(newunit=input_unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open trajectory file: " // trim(filename), ERR_FILE_IO)
        end if
        input_file_open = .true.
        
        ! Allocate module-level buffer
        if (allocated(values_buffer)) deallocate(values_buffer)
        allocate(values_buffer(1000))  ! Initial size, will grow if needed
        
        ! Check if user provided a header format in setup
        if (len_trim(user_header_line) > 0) then
            use_user_header = .true.
        end if
    end subroutine init_trajectory_reader

    !> @brief Get file extension from filename
    !>
    !> Extracts the extension from a filename (characters after the last dot).
    !>
    !> @param[in] filename Filename to analyze
    !> @return File extension (lowercase) or empty string if none found
    function get_file_extension(filename) result(ext)
        character(len=*), intent(in) :: filename
        character(len=10) :: ext
        integer :: i, dot_pos
        
        ext = ""
        dot_pos = 0
        
        ! Find the last dot in the filename
        do i = len_trim(filename), 1, -1
            if (filename(i:i) == '.') then
                dot_pos = i
                exit
            end if
        end do
        
        ! Extract extension if found
        if (dot_pos > 0) then
            ext = filename(dot_pos+1:)
        end if
    end function get_file_extension

    !> @brief Count frames in trajectory file
    !>
    !> Counts the number of frames in the trajectory file and optionally
    !> builds a frame offset table for fast seeking.
    !>
    !> @param[in] filename Path to trajectory file
    !> @return Number of frames in the file
    function count_trajectory_frames(filename) result(num_frames)
        character(len=*), intent(in) :: filename
        integer :: num_frames
        integer :: unit, io_stat, grep_stat, cmd_stat
        logical :: file_exists
        character(256) :: command, result_string
        character(len=32) :: temp_file
        character(len=32) :: grep_positions_file
        character(len=10) :: file_ext
        
        num_frames = 0
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) "Warning: Could not find file: ", trim(filename)
            return
        end if
        
        ! Determine file format based on extension
        file_ext = get_file_extension(filename)
        
        if (trim(file_ext) == "xyz") then
            ! For XYZ files, count occurrences of the atom count line
            
            ! Try using grep (faster on large files)
            temp_file = "frame_count_temp.txt"
            command = 'grep -c "^[0-9][0-9]*$" ' // trim(filename) // ' > ' // trim(temp_file)
            call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
            
            ! If grep worked, read the count
            if (grep_stat == 0 .and. cmd_stat == 0) then
                open(newunit=unit, file=temp_file, status='old', action='read', iostat=io_stat)
                if (io_stat == 0) then
                    read(unit, *, iostat=io_stat) num_frames
                    close(unit)
                end if
                
                ! Clean up temporary file
                call execute_command_line('rm ' // trim(temp_file), wait=.true.)
            else
                ! Fall back to manual counting
                num_frames = count_xyz_frames_direct(filename)
            end if
        else
            ! For LAMMPS files, use existing method
            
            temp_file = "frame_count_temp.txt"
            grep_positions_file = "frame_positions_temp.txt"
            
            ! Count frames and get byte offsets for each frame
            command = 'grep -b "ITEM: NUMBER OF ATOMS" ' // trim(filename) // ' > ' // trim(grep_positions_file)
            call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
            
            ! Get just the count
            command = 'grep -c "ITEM: NUMBER OF ATOMS" ' // trim(filename) // ' > ' // trim(temp_file)
            call execute_command_line(trim(command), exitstat=grep_stat, cmdstat=cmd_stat)
            
            ! If command execution fails or grep isn't available, fall back immediately
            if (grep_stat /= 0 .or. cmd_stat /= 0) then
                if (VERBOSE) write(*,*) "Using direct frame counting method"
                num_frames = count_lammps_frames_direct(filename)
                
                ! Clean up temporary files if they were created
                inquire(file=temp_file, exist=file_exists)
                if (file_exists) then
                    call execute_command_line('rm ' // trim(temp_file), wait=.true.)
                end if
                
                inquire(file=grep_positions_file, exist=file_exists)
                if (file_exists) then
                    call execute_command_line('rm ' // trim(grep_positions_file), wait=.true.)
                end if
                
                return
            end if
            
            ! Read the frame count
            open(newunit=unit, file=temp_file, status='old', action='read', iostat=io_stat)
            if (io_stat == 0) then
                read(unit, *, iostat=io_stat) num_frames
                close(unit)
            end if
            
            ! Parse the frame offsets if we have a valid count
            if (num_frames > 0) then
                call parse_frame_offsets(grep_positions_file, num_frames)
            end if
            
            ! Clean up temporary files
            call execute_command_line('rm ' // trim(temp_file), wait=.true.)
            call execute_command_line('rm ' // trim(grep_positions_file), wait=.true.)
        end if
    end function count_trajectory_frames
    
    !> @brief Count frames in XYZ file directly
    !>
    !> Falls back to direct file reading to count frames in an XYZ file
    !> when grep is not available or fails.
    !>
    !> @param[in] fname Path to XYZ file
    !> @return Number of frames in the file
    function count_xyz_frames_direct(fname) result(count)
        character(len=*), intent(in) :: fname
        integer :: count
        integer :: cnt_unit, cnt_stat, i
        character(256) :: cnt_line
        integer :: line_count, n_atoms
        
        count = 0
        open(newunit=cnt_unit, file=fname, status='old', action='read', iostat=cnt_stat)
        if (cnt_stat /= 0) return
        
        line_count = 0
        do
            read(cnt_unit, '(A)', iostat=cnt_stat) cnt_line
            if (cnt_stat /= 0) exit
            
            line_count = line_count + 1
            if (line_count == 1) then
                ! First line of frame contains atom count
                read(cnt_line, *, iostat=cnt_stat) n_atoms
                if (cnt_stat == 0) then
                    count = count + 1
                    ! Skip n_atoms + 1 lines (comment line + atom data)
                    do i = 1, n_atoms + 1
                        read(cnt_unit, '(A)', iostat=cnt_stat) cnt_line
                        if (cnt_stat /= 0) exit
                    end do
                    if (cnt_stat /= 0) exit
                    line_count = 0
                end if
            end if
        end do
        
        close(cnt_unit)
    end function count_xyz_frames_direct
    
    !> @brief Count frames in LAMMPS file directly
    !>
    !> Falls back to direct file reading to count frames in a LAMMPS file
    !> when grep is not available or fails.
    !>
    !> @param[in] fname Path to LAMMPS file
    !> @return Number of frames in the file
    function count_lammps_frames_direct(fname) result(count)
        character(len=*), intent(in) :: fname
        integer :: count
        integer :: cnt_unit, cnt_stat
        character(256) :: cnt_line
        
        count = 0
        open(newunit=cnt_unit, file=fname, status='old', action='read', iostat=cnt_stat)
        if (cnt_stat /= 0) return
        
        do
            read(cnt_unit, '(A)', iostat=cnt_stat) cnt_line
            if (cnt_stat /= 0) exit
            if (index(cnt_line, 'ITEM: NUMBER OF ATOMS') > 0) count = count + 1
        end do
        
        close(cnt_unit)
    end function count_lammps_frames_direct
    
    !> @brief Parse frame offsets from grep output
    !>
    !> Processes the output of 'grep -b' to build a table of file offsets
    !> for each frame, enabling fast seeking to specific frames.
    !>
    !> @param[in] filename File containing grep output with byte offsets
    !> @param[in] num_frames Number of frames to parse offsets for
    subroutine parse_frame_offsets(filename, num_frames)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: num_frames
        integer :: unit, io_stat, i, alloc_stat
        character(256) :: line
        integer(kind=8) :: offset
        
        ! Allocate array for frame offsets
        if (allocated(frame_offsets)) deallocate(frame_offsets)
        allocate(frame_offsets(num_frames), stat=alloc_stat)
        if (alloc_stat /= 0) then
            if (VERBOSE) write(*,*) "Warning: Could not allocate frame offsets array"
            have_frame_offsets = .false.
            return
        end if
        
        ! Open the grep output file
        open(newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            have_frame_offsets = .false.
            return
        end if
        
        ! Read each line and extract the byte offset
        do i = 1, num_frames
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Extract offset (format is "OFFSET:ITEM: NUMBER OF ATOMS")
            read(line(:index(line, ':')-1), *, iostat=io_stat) offset
            if (io_stat == 0) then
                frame_offsets(i) = offset
            else
                if (VERBOSE) write(*,*) "Warning: Error parsing frame offset at line", i
            end if
        end do
        
        close(unit)
        
        ! Check if we got all offsets
        if (i > num_frames) then
            have_frame_offsets = .true.
        else
            if (VERBOSE) write(*,*) "Warning: Only found", i-1, "frame offsets of", num_frames
            have_frame_offsets = (i-1 >= num_frames)
        end if
    end subroutine parse_frame_offsets

    !> @brief Seek to a specific frame in the trajectory file
    !>
    !> Positions the file pointer at the beginning of the specified frame.
    !> Uses direct seeking with frame offsets if available, otherwise
    !> seeks sequentially by reading through the file.
    !>
    !> @param[in] frame_number Frame to seek to (0-based, 0 means rewind)
    subroutine seek_to_frame(frame_number)
        integer, intent(in) :: frame_number
        integer :: i, j, io_stat, n_atoms_to_skip
        character(256) :: dummy_line
        logical :: success
        
        if (.not. input_file_open) then
            call handle_error("File not open for seeking", ERR_FILE_IO)
            return
        end if
        
        ! If frame number is 0 or negative, just rewind
        if (frame_number <= 0) then
            rewind(input_unit)
            return
        end if
        
        ! If we have frame offsets, use direct seeking for LAMMPS format
        if (file_format == FORMAT_LAMMPS .and. have_frame_offsets .and. &
            frame_number > 0 .and. frame_number <= size(frame_offsets)) then
            
            ! Seek to the position of the frame
            rewind(input_unit)
            
            ! More robust approach to positioning
            success = .false.
            
            ! Try direct positioning first
            read(input_unit, '(A)', advance='no', pos=frame_offsets(frame_number), iostat=io_stat) dummy_line
            if (io_stat == 0) then
                ! Check we're at the right frame marker
                if (index(dummy_line, 'ITEM: NUMBER OF ATOMS') > 0) then
                    success = .true.
                end if
            end if
            
            ! If direct positioning failed or gave wrong position, fall back to sequential
            if (.not. success) then
                if (VERBOSE) write(*,*) "Direct frame seeking failed, using sequential method"
                rewind(input_unit)
                
                ! Skip frames sequentially
                do i = 1, frame_number - 1
                    ! Find the start of a frame
                    do
                        read(input_unit, '(A)', iostat=io_stat) dummy_line
                        if (io_stat /= 0) then
                            call handle_error("Reached end of file while seeking frame " // &
                                           trim(integer_to_string(frame_number)), ERR_FILE_IO)
                            return
                        end if
                        
                        if (index(dummy_line, 'ITEM: NUMBER OF ATOMS') > 0) exit
                    end do
                    
                    ! Read atom count
                    read(input_unit, *, iostat=io_stat) n_atoms_to_skip
                    if (io_stat /= 0) then
                        call handle_error("Error reading atom count while seeking", ERR_FILE_IO)
                        return
                    end if
                    
                    ! Find atom section
                    do
                        read(input_unit, '(A)', iostat=io_stat) dummy_line
                        if (io_stat /= 0) then
                            call handle_error("Reached end of file while seeking", ERR_FILE_IO)
                            return
                        end if
                        
                        if (index(dummy_line, 'ITEM: ATOMS') > 0) exit
                    end do
                    
                    ! Skip atom lines
                    do j = 1, n_atoms_to_skip
                        read(input_unit, '(A)', iostat=io_stat) dummy_line
                        if (io_stat /= 0) then
                            call handle_error("Reached end of file while seeking", ERR_FILE_IO)
                            return
                        end if
                    end do
                end do
            end if
            
            return
        end if
        
        ! For XYZ format or if no offsets, fall back to sequential reading
        rewind(input_unit)
        
        if (file_format == FORMAT_LAMMPS) then
            ! Sequential seek for LAMMPS format
            do i = 1, frame_number - 1
                ! Find the start of a frame
                do
                    read(input_unit, '(A)', iostat=io_stat) dummy_line
                    if (io_stat /= 0) then
                        call handle_error("Reached end of file while seeking frame " // &
                                        trim(integer_to_string(frame_number)), ERR_FILE_IO)
                        return
                    end if
                    
                    if (index(dummy_line, 'ITEM: NUMBER OF ATOMS') > 0) exit
                end do
                
                ! Read atom count
                read(input_unit, *, iostat=io_stat) n_atoms_to_skip
                if (io_stat /= 0) then
                    call handle_error("Error reading atom count while seeking", ERR_FILE_IO)
                    return
                end if
                
                ! Find atom section
                do
                    read(input_unit, '(A)', iostat=io_stat) dummy_line
                    if (io_stat /= 0) then
                        call handle_error("Reached end of file while seeking", ERR_FILE_IO)
                        return
                    end if
                    
                    if (index(dummy_line, 'ITEM: ATOMS') > 0) exit
                end do
                
                ! Skip atom lines
                do j = 1, n_atoms_to_skip
                    read(input_unit, '(A)', iostat=io_stat) dummy_line
                    if (io_stat /= 0) then
                        call handle_error("Reached end of file while seeking", ERR_FILE_IO)
                        return
                    end if
                end do
            end do
        else
            ! For XYZ format
            do i = 1, frame_number - 1
                ! Read atom count
                read(input_unit, '(A)', iostat=io_stat) dummy_line
                if (io_stat /= 0) then
                    call handle_error("Reached end of file before frame " // trim(integer_to_string(frame_number)), &
                                    ERR_FILE_IO)
                    return
                end if
                
                ! Convert atom count to integer
                read(dummy_line, *, iostat=io_stat) n_atoms_to_skip
                if (io_stat /= 0) then
                    call handle_error("Error reading atom count in XYZ file", ERR_FILE_IO)
                    return
                end if
                
                ! Skip comment line
                read(input_unit, '(A)', iostat=io_stat) dummy_line
                if (io_stat /= 0) exit
                
                ! Skip atom data lines
                do j = 1, n_atoms_to_skip
                    read(input_unit, '(A)', iostat=io_stat) dummy_line
                    if (io_stat /= 0) exit
                end do
                if (io_stat /= 0) then
                    call handle_error("Reached end of file while skipping frame", ERR_FILE_IO)
                    return
                end if
            end do
        end if
    end subroutine seek_to_frame
    
    !> @brief Convert integer to string
    !>
    !> Helper function to convert an integer to a string.
    !>
    !> @param[in] int_value Integer value to convert
    !> @return String representation of the integer
    function integer_to_string(int_value) result(str)
        integer, intent(in) :: int_value
        character(len=20) :: str
        
        write(str, '(I0)') int_value
    end function integer_to_string

    !> @brief Read a frame from the trajectory file
    !>
    !> Reads a single frame from the trajectory file into the provided arrays.
    !> Automatically dispatches to the appropriate format-specific reader.
    !>
    !> @param[out] coords Atom coordinates
    !> @param[out] atom_types Atom type indices
    !> @param[out] elements Atom element names
    !> @param[out] box_length Box dimensions in x, y, z
    !> @param[inout] frame_number Current frame number
    !> @param[out] eof Flag indicating end-of-file
    !> @param[in] atom_info Array of atom type information
    subroutine read_trajectory_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        real(dp), intent(out) :: coords(:,:)
        integer, intent(out) :: atom_types(:)
        character(len=*), intent(out) :: elements(:)
        real(dp), intent(out) :: box_length(3)
        integer, intent(inout) :: frame_number
        logical, intent(out) :: eof
        type(atom_type_info), intent(in) :: atom_info(:)
        
        if (file_format == FORMAT_LAMMPS) then
            call read_lammps_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        else
            call read_xyz_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        end if
    end subroutine read_trajectory_frame

    !> @brief Read a frame from a LAMMPS trajectory file
    !>
    !> Specialized reader for the LAMMPS dump file format.
    !>
    !> @param[out] coords Atom coordinates
    !> @param[out] atom_types Atom type indices
    !> @param[out] elements Atom element names
    !> @param[out] box_length Box dimensions in x, y, z
    !> @param[inout] frame_number Current frame number
    !> @param[out] eof Flag indicating end-of-file
    !> @param[in] atom_info Array of atom type information
    subroutine read_lammps_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        real(dp), intent(out) :: coords(:,:)
        integer, intent(out) :: atom_types(:)
        character(len=*), intent(out) :: elements(:)
        real(dp), intent(out) :: box_length(3)
        integer, intent(inout) :: frame_number
        logical, intent(out) :: eof
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: n_atoms_frame, atom_id, io_stat, i, timestep
        real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
        character(512) :: line
        
        eof = .false.
        
        ! Look for frame start
        do
            read(input_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                eof = .true.
                return
            end if
            
            if (index(line, 'ITEM: NUMBER OF ATOMS') > 0) then
                read(input_unit, *) n_atoms_frame
                if (n_atoms_frame /= size(coords,1)) then
                    write(*,*) "Warning: Atom count mismatch -", &
                              "Expected:", size(coords,1), &
                              "Found:", n_atoms_frame
                end if
                exit
            end if
        end do
        
        ! Read box bounds
        read(input_unit, '(A)') line  ! ITEM: BOX BOUNDS
        read(input_unit, *) xlo, xhi
        read(input_unit, *) ylo, yhi
        read(input_unit, *) zlo, zhi
        
        ! Calculate box lengths
        box_length(1) = xhi - xlo
        box_length(2) = yhi - ylo
        box_length(3) = zhi - zlo
        
        ! Read atoms header and setup column mapping
        read(input_unit, '(A)') line
        
        ! If user provided a header line, use it instead of the one in the file
        if (use_user_header) then
            line = user_header_line
        end if
        
        ! Setup column mapping for first frame or if using user header
        if (frame_number == 0 .or. use_user_header) then
            call setup_column_mapping(line)
            header_shown = .true.  ! Mark that we've shown the header info
        end if
        
        ! Ensure values_buffer is large enough
        if (.not. allocated(values_buffer) .or. size(values_buffer) < col_map%max_cols) then
            if (allocated(values_buffer)) deallocate(values_buffer)
            allocate(values_buffer(col_map%max_cols))
        end if
        
        ! Read atom data
        do i = 1, n_atoms_frame
            read(input_unit, *, iostat=io_stat) values_buffer(1:col_map%max_cols)
            if (io_stat /= 0) then
                write(*,*) "Error reading atom", i, "in frame", frame_number
                call handle_error("Error reading coordinates", ERR_FILE_IO)
            end if
            
            atom_id = nint(values_buffer(col_map%id))
            if (atom_id < 1 .or. atom_id > size(coords,1)) then
                write(*,*) "Invalid atom ID:", atom_id, "at line", i
                call handle_error("Invalid atom ID in trajectory", ERR_INCONSISTENT_DATA)
            end if
            
            coords(atom_id,1) = values_buffer(col_map%x)
            coords(atom_id,2) = values_buffer(col_map%y)
            coords(atom_id,3) = values_buffer(col_map%z)
            atom_types(atom_id) = nint(values_buffer(col_map%type))
            
            ! Set element name from atom_info instead of generic "Type-X"
            if (atom_types(atom_id) >= 1 .and. atom_types(atom_id) <= size(atom_info)) then
                elements(atom_id) = atom_info(atom_types(atom_id))%name
            else
                write(elements(atom_id),'(A,I0)') 'Type', atom_types(atom_id)
            end if
        end do
        
        frame_number = frame_number + 1
    end subroutine read_lammps_frame

    !> @brief Read a frame from an XYZ trajectory file
    !>
    !> Specialized reader for the XYZ file format.
    !>
    !> @param[out] coords Atom coordinates
    !> @param[out] atom_types Atom type indices
    !> @param[out] elements Atom element names
    !> @param[out] box_length Box dimensions in x, y, z
    !> @param[inout] frame_number Current frame number
    !> @param[out] eof Flag indicating end-of-file
    !> @param[in] atom_info Array of atom type information
    subroutine read_xyz_frame(coords, atom_types, elements, box_length, frame_number, eof, atom_info)
        real(dp), intent(out) :: coords(:,:)
        integer, intent(out) :: atom_types(:)
        character(len=*), intent(out) :: elements(:)
        real(dp), intent(out) :: box_length(3)
        integer, intent(inout) :: frame_number
        logical, intent(out) :: eof
        type(atom_type_info), intent(in) :: atom_info(:)
        
        integer :: n_atoms_frame, i, io_stat
        character(256) :: line, element_str
        real(dp) :: x, y, z
        logical :: box_info_found
        
        eof = .false.
        box_info_found = .false.
        
        ! Read atom count for this frame
        read(input_unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) then
            eof = .true.
            return
        end if
        
        ! Extract atom count
        read(line, *, iostat=io_stat) n_atoms_frame
        if (io_stat /= 0) then
            call handle_error("Error reading atom count in XYZ file", ERR_FILE_IO)
            return
        end if
        
        if (n_atoms_frame /= size(coords,1)) then
            write(*,*) "Warning: Atom count mismatch -", &
                      "Expected:", size(coords,1), &
                      "Found:", n_atoms_frame
        end if
        
        ! Read comment line, which might contain box information
        read(input_unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) then
            eof = .true.
            return
        end if
        
        ! Try to extract box information from the comment line
        ! Format: "Box: x_length y_length z_length" or similar
        if (index(line, 'Box') > 0 .or. index(line, 'box') > 0) then
            ! Extract box dimensions from the comment line
            read(line(index(line, ':')+1:), *, iostat=io_stat) box_length(1), box_length(2), box_length(3)
            if (io_stat == 0) then
                box_info_found = .true.
            end if
        end if
        
        ! If no box info found, use previous box_length or defaults
        if (.not. box_info_found) then
            ! Keep the existing box_length values (already set from data file or previous frame)
            if (VERBOSE .and. frame_number == 0) then
                write(*,*) "No box information found in XYZ comment line, using values from data file."
            end if
        end if
        
        ! Read atom data
        do i = 1, n_atoms_frame
            read(input_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                write(*,*) "Error reading atom", i, "in frame", frame_number
                call handle_error("Error reading coordinates in XYZ file", ERR_FILE_IO)
                return
            end if
            
            ! Parse line: first field is element, next three are x, y, z
            read(line, *, iostat=io_stat) element_str, x, y, z
            if (io_stat /= 0) then
                write(*,*) "Error parsing atom data in XYZ file, line:", trim(line)
                call handle_error("Error parsing atom data in XYZ file", ERR_FILE_IO)
                return
            end if
            
            ! Store coordinates
            coords(i,1) = x
            coords(i,2) = y
            coords(i,3) = z
            
            ! Determine atom type from element name
            elements(i) = trim(element_str)
            atom_types(i) = get_atom_type_from_element(element_str, atom_info)
        end do
        
        frame_number = frame_number + 1
    end subroutine read_xyz_frame
    
    !> @brief Map element name to atom type
    !>
    !> Determines the atom type index based on element name, using
    !> the atom_info array to match names to type IDs.
    !>
    !> @param[in] element Element name to look up
    !> @param[in] atom_info Array of atom type information
    !> @return Atom type index
    function get_atom_type_from_element(element, atom_info) result(type_id)
        character(len=*), intent(in) :: element
        type(atom_type_info), intent(in) :: atom_info(:)
        integer :: type_id, i, io_stat
        
        ! Try to find a matching element name in atom_info
        type_id = 0
        do i = 1, size(atom_info)
            if (trim(element) == trim(atom_info(i)%name)) then
                type_id = i
                return
            end if
        end do
        
        ! If no match found, assume it's an unknown element
        type_id = 1  ! Default to first type
        
        ! Try to extract a number if the element is in "TypeX" format
        if (element(1:4) == 'Type' .or. element(1:4) == 'type') then
            read(element(5:), *, iostat=io_stat) i
            if (io_stat == 0 .and. i >= 1 .and. i <= size(atom_info)) then
                type_id = i
            end if
        end if
    end function get_atom_type_from_element

    !> @brief Set up column mapping for trajectory data
    !>
    !> Parses the header line of a LAMMPS dump file to determine which
    !> columns contain atom IDs, types, and coordinates.
    !>
    !> @param[in] header_line Header line to parse
    subroutine setup_column_mapping(header_line)
        character(len=*), intent(in) :: header_line
        integer :: i, pos, prev_pos, n_cols
        character(len=32) :: column_name
        
        ! Reset column mapping
        col_map = column_map_type()
        
        ! Skip "ITEM: ATOMS" prefix if present
        if (index(header_line, "ITEM: ATOMS") > 0) then
            prev_pos = index(header_line, "ITEM: ATOMS") + 11
        else
            prev_pos = 1
        end if
        
        ! Enhanced parsing to find column indices
        n_cols = 0
        pos = prev_pos
        
        ! First count columns
        do while (pos <= len_trim(header_line))
            ! Skip spaces
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) == ' ')
                pos = pos + 1
            end do
            
            ! If we're at the end, break
            if (pos > len_trim(header_line)) exit
            
            ! Found a non-space, this is a column
            n_cols = n_cols + 1
            
            ! Skip the column name
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) /= ' ')
                pos = pos + 1
            end do
        end do
        
        if (VERBOSE .and. .not. shown_column_info) then
            write(*,*) "Found", n_cols, "columns in header line"
        end if
        
        ! Now map column indices
        pos = prev_pos
        do i = 1, n_cols
            ! Skip spaces
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) == ' ')
                pos = pos + 1
            end do
            
            ! If we're at the end, break
            if (pos > len_trim(header_line)) exit
            
            ! Extract column name
            prev_pos = pos
            do while (pos <= len_trim(header_line) .and. header_line(pos:pos) /= ' ')
                pos = pos + 1
            end do
            
            column_name = adjustl(header_line(prev_pos:pos-1))
            
            ! Map column - more robust checking for variations
            if (column_name == 'id' .or. column_name == 'Id' .or. column_name == 'ID') then
                col_map%id = i
            else if (column_name == 'type' .or. column_name == 'Type' .or. column_name == 'TYPE') then
                col_map%type = i
            else if (column_name == 'x' .or. column_name == 'X' .or. column_name == 'xu') then
                col_map%x = i
            else if (column_name == 'y' .or. column_name == 'Y' .or. column_name == 'yu') then
                col_map%y = i
            else if (column_name == 'z' .or. column_name == 'Z' .or. column_name == 'zu') then
                col_map%z = i
            end if
        end do
        
        ! Set maximum columns needed
        col_map%max_cols = max(col_map%id, col_map%type, col_map%x, col_map%y, col_map%z)
        
        ! Verify required columns
        if (col_map%id < 0 .or. col_map%type < 0 .or. &
            col_map%x < 0 .or. col_map%y < 0 .or. col_map%z < 0) then
            call handle_error("Missing required columns in trajectory file", ERR_FILE_IO)
        end if
    end subroutine setup_column_mapping

    !> @brief Clean up trajectory reader resources
    !>
    !> Closes the trajectory file and deallocates module arrays.
    subroutine cleanup_trajectory_reader()
        if (input_file_open) then
            close(input_unit)
            input_file_open = .false.
        end if
        
        if (allocated(frame_offsets)) deallocate(frame_offsets)
        if (allocated(values_buffer)) deallocate(values_buffer)
        have_frame_offsets = .false.
        
        ! Reset flags for next use
        header_shown = .false.
        shown_column_info = .false.
    end subroutine cleanup_trajectory_reader

    !> @brief Set custom header line for trajectory parsing
    !>
    !> Allows specifying a custom header format for parsing trajectory files,
    !> useful when the file format is non-standard or needs customization.
    !>
    !> @param[in] header Custom header line to use for column mapping
    subroutine set_user_header_line(header)
        character(len=*), intent(in) :: header
        
        if (len_trim(header) > 0) then
            user_header_line = header
            use_user_header = .true.
        end if
    end subroutine set_user_header_line

end module trajectory_io_mod
