!> @brief Main program for atomic coordination number analysis
program coordination_analysis
    use types_mod
    use config_mod
    use cell_list_mod
    use coordination_mod
    use io_mod
    use error_mod
    use atom_selection_mod
    use omp_lib
    use angle_mod
    implicit none

    !===============================================================================
    ! System variables
    !===============================================================================
    integer :: n_atoms = 0                   ! Number of atoms in the system
    integer :: n_types = 0                   ! Number of atom types specified in setup
    integer :: n_types_actual = 0            ! Number of atom types found in data file
    type(atom_type_info), allocatable :: atom_info(:)  ! Array of atom type information
    real(dp) :: box_length(3)                ! Box dimensions in x, y, z (Å)
    
    !===============================================================================
    ! Atomic data arrays
    !===============================================================================
    real(dp), allocatable :: coords(:,:)     ! Atom coordinates [atom_index, dimension]
    integer, allocatable :: atom_types(:)    ! Atom type indices
    character(len=8), allocatable :: elements(:)  ! Atom element names
    
    !===============================================================================
    ! Analysis variables
    !===============================================================================
    integer :: frame = 0                     ! Current frame number
    integer :: total_frames = 0              ! Total number of frames to process
    logical :: eof = .false.                 ! Flag indicating end of trajectory file
    real(dp) :: start_time                   ! Analysis start time
    real(dp) :: end_time                     ! Analysis end time
    logical, allocatable :: include_mask(:)  ! Mask indicating which atom types to include
    
    !===============================================================================
    ! Local variables
    !===============================================================================
    integer :: io_stat                       ! Status code for I/O operations
    integer :: num_threads                   ! Number of OpenMP threads to use
    character(len=32) :: arg                 ! Command-line argument buffer
    logical :: setup_file_read = .false.     ! Flag if setup file was successfully read
    logical :: atom_selection_initialized = .false.  ! Track if atom selection was initialized
    logical :: trajectory_header_set = .false.       ! Track if trajectory header was set
    real(dp) :: prev_box_length(3) = 0.0_dp  ! Track previous box dimensions
    logical :: box_changed = .false.         ! Flag to track if box dimensions changed
    
    ! Box dimension statistics variables
    integer :: box_change_count = 0           ! Counter for box dimension changes
    real(dp) :: min_box_volume = 0.0_dp       ! Minimum box volume
    real(dp) :: max_box_volume = 0.0_dp       ! Maximum box volume
    real(dp) :: total_box_volume = 0.0_dp     ! Accumulated box volume for average
    real(dp) :: first_box_volume = 0.0_dp     ! Initial box volume
    real(dp) :: last_box_volume = 0.0_dp      ! Final box volume
    real(dp) :: box_volume                    ! Current box volume
    
    ! Variable for type counts and percentages
    integer :: type_count, i, j, k, type1, type2, types_found, temp_id
    real(dp) :: percentage
    integer, allocatable :: type_counts(:)
    logical, allocatable :: explicit_pairs(:)   ! Track which pairs were explicitly defined
    
    ! Date and time variables
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer, dimension(8) :: values

    !===============================================================================
    ! Print program header with time
    !===============================================================================
    call date_and_time(date, time, zone, values)
    write(*,'(A)') " ==========================================================================="
    write(*,'(A,4X,A4,A1,A2,A1,A2,1X,A2,A1,A2,A1,A2)') &
              " COORDINATION ANALYSIS   ", date(1:4), "-", date(5:6), "-", date(7:8), &
              time(1:2), ":", time(3:4), ":", time(5:6)
    write(*,'(A)') " ==========================================================================="
    write(*,*)

    !===============================================================================
    ! Initialize configuration
    !===============================================================================
    call initialize_config()
    call read_config()
    
    !===============================================================================
    ! Read system dimensions from data file
    !===============================================================================
    if (.not. use_trajectory_box) then
        call read_data_file_dims_only(DATA_FILE, n_atoms, n_types_actual, box_length)
        
        write(*,'(A)') " SYSTEM CONFIGURATION"
        write(*,'(A)') " -------------------"
        write(*,'(A,I0)') " Found atoms:                ", n_atoms 
        write(*,'(A,I0)') " Found atom types:           ", n_types_actual
        write(*,*)
        write(*,'(A)') " System dimensions:"
        write(*,'(A,F10.4)') "   X length:                 ", box_length(1)
        write(*,'(A,F10.4)') "   Y length:                 ", box_length(2)
        write(*,'(A,F10.4)') "   Z length:                 ", box_length(3)
        write(*,*)
        
        ! Calculate box volume
        box_volume = box_length(1) * box_length(2) * box_length(3)
        write(*,'(A,F10.1,A)') " Box volume:                 ", box_volume, " cubic Angstroms"
        write(*,*)
    else
        ! When using trajectory box, we still need to get atom count and types
        call read_data_file_dims_only(DATA_FILE, n_atoms, n_types_actual, box_length)
        
        write(*,'(A)') " SYSTEM CONFIGURATION"
        write(*,'(A)') " -------------------"
        write(*,'(A,I0)') " Found atoms:                ", n_atoms 
        write(*,'(A,I0)') " Found atom types:           ", n_types_actual
        write(*,*)
        write(*,'(A)') " System dimensions will be read from trajectory file"
        write(*,*)
    end if
    
    !===============================================================================
    ! Read setup file
    !===============================================================================
    setup_file_read = read_setup_file(atom_info, n_types, pairs, n_pairs)
    
    if (setup_file_read) then
        ! Configure OpenMP based on setup
        if (requested_cores > 0) then
            ! Use number of cores from setup file
            num_threads = min(requested_cores, omp_get_max_threads())
        else
            ! Use all available cores
            num_threads = omp_get_max_threads()
        end if
        
        ! Set OpenMP parameters
        call omp_set_num_threads(num_threads)
        call omp_set_schedule(omp_sched_guided, 0)
    else
        ! If no setup file, use all cores
        num_threads = omp_get_max_threads()
        call omp_set_num_threads(num_threads)
        call omp_set_schedule(omp_sched_guided, 0)
    end if
    
    !===============================================================================
    ! Allocate system arrays
    !===============================================================================
    allocate(coords(n_atoms,3), atom_types(n_atoms), elements(n_atoms), &
             include_mask(n_types), stat=io_stat)
    if (io_stat /= 0) call handle_error("Failed to allocate arrays", ERR_ALLOCATION)
    
    ! Initialize arrays
    coords = 0.0_dp
    atom_types = 0
    elements = ''
    include_mask = atom_info%include
    
    ! Allocate and initialize type counts and explicit pairs arrays
    allocate(type_counts(n_types), explicit_pairs(n_pairs))
    type_counts = 0
    explicit_pairs = .false.
    
    ! Mark pairs that were explicitly defined in setup.txt
    ! We'll consider a pair as explicitly defined if it doesn't use the default cutoff
    ! or if a self-interaction pair (same type) which is always defined
    do i = 1, n_pairs
        if (pairs(i)%cutoff /= DEFAULT_CUTOFF .or. pairs(i)%type1_id == pairs(i)%type2_id) then
            explicit_pairs(i) = .true.
        end if
    end do
    
    ! Read first frame to get atom types
    if (setup_file_read) then
        call read_trajectory_frame_noupdate(coords, atom_types, elements, box_length, &
                                      frame, eof, atom_info, input_trajectory)
    else
        call read_trajectory_frame_noupdate(coords, atom_types, elements, box_length, &
                                      frame, eof, atom_info, INPUT_FILE)
    end if
    
    ! If using trajectory box dimensions, use them from the first frame
    if (use_trajectory_box) then
        write(*,'(A)') " System dimensions from trajectory:"
        write(*,'(A,F10.4)') "   X length:                 ", box_length(1)
        write(*,'(A,F10.4)') "   Y length:                 ", box_length(2)
        write(*,'(A,F10.4)') "   Z length:                 ", box_length(3)
        write(*,*)
        
        ! Calculate box volume
        box_volume = box_length(1) * box_length(2) * box_length(3)
        first_box_volume = box_volume
        min_box_volume = box_volume
        max_box_volume = box_volume
        total_box_volume = box_volume
        
        write(*,'(A,F10.1,A)') " Box volume:                 ", box_volume, " cubic Angstroms"
        write(*,*)
    end if
    
    ! Count atoms by type
    do i = 1, n_atoms
        if (atom_types(i) >= 1 .and. atom_types(i) <= n_types) then
            type_counts(atom_types(i)) = type_counts(atom_types(i)) + 1
        end if
    end do
    
    ! Print atom composition
    write(*,'(A)') " Atom composition:"
    types_found = 0
    do i = 1, n_types
        ! Only show types that have atoms and valid names
        if (atom_info(i)%include .and. type_counts(i) > 0 .and. len_trim(atom_info(i)%name) > 0) then
            types_found = types_found + 1
            type_count = type_counts(i)
            percentage = real(type_count) / real(n_atoms) * 100.0
            
            ! Format with consistent spacing
            write(*,'(A,I0,A,A,A,I0,A,F6.2,A)') "   Type ", i, " (", &
                trim(atom_info(i)%name), "):             ", type_count, " atoms (", percentage, "%)"
        end if
    end do

    ! If no valid types were found, show a warning
    if (types_found == 0) then
        write(*,'(A)') "   Warning: No valid atom types found with non-zero counts"
    end if
    write(*,*)
    
    !===============================================================================
    ! Set up trajectory processing
    !===============================================================================
    write(*,'(A)') " ANALYSIS SETTINGS"
    write(*,'(A)') " ----------------"
    
    ! Print atom selection info
    if (len_trim(selected_atoms) > 0) then
        write(*,'(A,A)') " Atom selection:            ", trim(selected_atoms)
    else
        write(*,'(A)') " Atom selection:            all"
    end if
    
    ! Print frame selection info
    if (start_frame > 0) then
        write(*,'(A,I0,A,A)') " Frame selection:           ", start_frame, " to ", &
                             trim(end_frame_str)
    else
        write(*,'(A,A,A,A)') " Frame selection:           beginning to ", trim(end_frame_str)
    end if
    
    ! Print calculation method info
    if (use_cell_list) then
        write(*,'(A)') " Calculation method:        Cell List"
        write(*,'(A,I0)') " Cell update frequency:     ", cell_update_freq
        if (selective_rebuild) then
            write(*,'(A)') " Selective cell rebuilding: True"
        else
            write(*,'(A)') " Selective cell rebuilding: False"
        end if
    else
        write(*,'(A)') " Calculation method:        Direct (O(n²))"
        write(*,'(A)') " Note: Direct method may be slower but useful for validation"
    end if
    
    if (use_trajectory_box) then
        write(*,'(A)') " Box dimensions:           From trajectory"
    else
        write(*,'(A)') " Box dimensions:           From data file"
    end if
    
    write(*,'(A,I0)') " OpenMP threads:            ", num_threads
    write(*,*)
    
    write(*,'(A)') " Trajectory format:         LAMMPS"
    if (len_trim(trajectory_header_format) > 0) then
        write(*,'(A,A)') " Trajectory header:         ", trim(trajectory_header_format)
    end if
    
    ! Set custom header format if provided
    if (len_trim(trajectory_header_format) > 0 .and. .not. trajectory_header_set) then
        call set_trajectory_header(trajectory_header_format)
        trajectory_header_set = .true.
    end if
    
    ! Initialize atom selection if specified
    if (len_trim(selected_atoms) > 0 .and. .not. atom_selection_initialized) then
        call initialize_atom_selection(n_atoms)
        call parse_atom_selection(selected_atoms, n_atoms)
        atom_selection_initialized = .true.
    end if
    
    ! Add column mapping section
    write(*,'(A)') " Column mapping:"
    write(*,'(A)') "   ID column:               1"
    write(*,'(A)') "   Type column:             2"
    write(*,'(A)') "   X column:                3"
    write(*,'(A)') "   Y column:                4"
    write(*,'(A)') "   Z column:                5"
    write(*,*)
    
    !===============================================================================
    ! Count frames in trajectory file
    !===============================================================================
    if (setup_file_read) then
        total_frames = count_trajectory_frames(input_trajectory)
    else 
        total_frames = count_trajectory_frames(INPUT_FILE)
    end if
    
    if (setup_file_read) then
        if (end_frame > 0) then
            total_frames = min(total_frames, end_frame)
        end if
    end if
    
    !===============================================================================
    ! Main program execution - normal analysis
    !===============================================================================
    write(*,'(A)') " PAIR CUTOFFS"
    write(*,'(A)') " -----------"
    
    ! Display only pair cutoffs explicitly defined in setup file
    k = 0 ! Counter for number of explicit pairs
    do i = 1, n_pairs
        ! Extract types for this pair
        type1 = pairs(i)%type1_id
        type2 = pairs(i)%type2_id
        
        ! Skip pairs where either type has 0 atoms or empty names
        if (type_counts(type1) == 0 .or. type_counts(type2) == 0) cycle
        if (len_trim(atom_info(type1)%name) == 0 .or. len_trim(atom_info(type2)%name) == 0) cycle
        
        ! Skip pairs that weren't explicitly defined in setup.txt
        if (.not. explicit_pairs(i)) cycle
        
        ! Use a safer consistent format that works for all cases
        write(*,'(A,I0,A,A,A,I0,A,A,A,F5.3,A)') &
            " Type ", type1, " (", trim(atom_info(type1)%name), ") - Type ", &
            type2, " (", trim(atom_info(type2)%name), "):   ", &
            pairs(i)%cutoff, " Angstrom"
        
        k = k + 1
    end do
    
    ! Display message if no explicit pairs were found
    if (k == 0) then
        write(*,'(A)') " No explicit pair cutoffs defined in setup file"
    end if
    write(*,*)
    
    ! Set up analysis - Initialize coordination
    call initialize_coordination(n_atoms, n_types, atom_info)
    
    ! Initialize cell list only if using that method
    if (use_cell_list) then
        call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
        
        ! Display cell list parameters section
        write(*,'(A)') " CELL LIST PARAMETERS"
        write(*,'(A)') " ------------------"
        write(*,'(A,F10.4,A)') " Cell size:                 ", cell_size, " Angstrom"
        write(*,'(A,I0,A,I0,A,I0)') " Grid dimensions:           ", n_cells_x, " x ", n_cells_y, " x ", n_cells_z
        write(*,'(A,I0)') " Total cells:               ", n_cells_x * n_cells_y * n_cells_z
        write(*,*)
    end if
    
    ! Initialize angle calculation if enabled
    if (calculate_angles) then
        call init_angle_analysis(n_atoms, n_types, atom_info)
    end if
    
    ! Initialize I/O
    if (setup_file_read) then
        call initialize_io(input_trajectory, output_data_file, total_frames, start_frame)
        frame = start_frame
    else
        call initialize_io(INPUT_FILE, OUTPUT_FILE, total_frames)
    end if
    
    ! Display trajectory analysis info
    write(*,'(A)') " TRAJECTORY ANALYSIS"
    write(*,'(A)') " ------------------"
    write(*,'(A,I0)') " Frames detected:           ", total_frames
    write(*,'(A,I0)') " Frames to process:         ", total_frames
    
    ! Display calculation method
    if (use_cell_list) then
        write(*,'(A)') " Method:                    Cell-List"
    else
        write(*,'(A)') " Method:                    Direct O(n²)"
    end if
    write(*,*)
    
    write(*,'(A)') " CALCULATION PROGRESS"
    write(*,'(A)') " -------------------"
    
    !===========================================================================
    ! Process trajectory frames
    !===========================================================================
    frame = 0  ! Reset frame counter
    start_time = omp_get_wtime()

    ! Initialize previous box dimensions
    prev_box_length = box_length

    do while (.not. eof)
        call read_trajectory_frame(coords, atom_types, elements, box_length, &
                                 frame, eof, atom_info)
        if (eof) exit
        
        ! Check if box dimensions have changed when using trajectory box
        box_changed = .false.
        if (use_trajectory_box) then
            if (abs(box_length(1) - prev_box_length(1)) > 1.0e-6_dp .or. &
                abs(box_length(2) - prev_box_length(2)) > 1.0e-6_dp .or. &
                abs(box_length(3) - prev_box_length(3)) > 1.0e-6_dp) then
                box_changed = .true.
                prev_box_length = box_length
                
                ! Update box statistics
                box_change_count = box_change_count + 1
                
                ! Calculate current box volume
                box_volume = box_length(1) * box_length(2) * box_length(3)
                last_box_volume = box_volume  ! Will contain the final volume at end
                
                ! Update min/max statistics
                min_box_volume = min(min_box_volume, box_volume)
                max_box_volume = max(max_box_volume, box_volume)
                
                total_box_volume = total_box_volume + box_volume
            end if
        end if
        
        ! Calculate coordination, only forcing rebuild if box actually changed
        ! Only calculate coordination for explicitly defined pairs
        call calculate_coordination(coords, atom_types, n_atoms, box_length, &
                                  frame, include_mask, box_changed)
        
        ! Calculate angles if enabled - with quiet mode for all frames except first
        if (calculate_angles) then
            if (frame == 0) then
                ! First frame shows full output
                call compute_angles(coords, atom_types, n_atoms, box_length, frame, atom_info)
            else
                ! Subsequent frames use quiet mode
                call compute_angles(coords, atom_types, n_atoms, box_length, frame, atom_info, .true.)
            end if
        end if
        
        call write_output_frame(frame, elements, atom_types, coord_numbers, n_atoms, atom_info)
        call update_progress(frame, total_frames)
        
        if (setup_file_read .and. end_frame > 0 .and. frame >= end_frame) exit
    end do
    
    end_time = omp_get_wtime()
    
    !===========================================================================
    ! Generate final report and cleanup
    !===========================================================================
    write(*,*)
    write(*,'(A)') " ==========================================================================="
    write(*,'(A)') "                            ANALYSIS COMPLETE"
    write(*,'(A)') " ==========================================================================="
    write(*,*)
    
    write(*,'(A)') " CALCULATION PERFORMANCE"
    write(*,'(A)') " ----------------------"
    write(*,'(A,I0)') " Processed frames:          ", frame
    write(*,'(A,I0,A)') " Processed atoms:           ", n_atoms, " per frame"
    write(*,'(A,F10.3,A)') " Total time:                ", end_time - start_time, " seconds"
    
    ! Calculate frames per second
    write(*,'(A,F10.3,A)') " Performance:               ", &
                          real(frame)/max(end_time - start_time, 1.0e-6_dp), " frames/second"
    
    ! Display calculation method
    if (use_cell_list) then
        write(*,'(A)') " Method used:               Cell-List"
    else
        write(*,'(A)') " Method used:               Direct O(n²)"
    end if
    write(*,*)
    
    ! Display box dimension statistics if trajectory box was used
    if (use_trajectory_box) then
        write(*,'(A)') " BOX DIMENSION STATISTICS"
        write(*,'(A)') " ----------------------"
        write(*,'(A,I0,A,I0,A)') " Box dimension changes: ", box_change_count, " (", &
                              nint(real(box_change_count)/max(1,frame)*100.0), "% of frames)"
        
        if (frame > 0) then
            write(*,'(A,F10.3,A)') " Initial box volume: ", first_box_volume, " cubic Angstroms"
            write(*,'(A,F10.3,A)') " Final box volume:   ", last_box_volume, " cubic Angstroms"
            write(*,'(A,F10.3,A)') " Change percentage:  ", &
                (last_box_volume-first_box_volume)/first_box_volume*100.0, "%"
            write(*,'(A,F10.3,A)') " Minimum box volume: ", min_box_volume, " cubic Angstroms"
            write(*,'(A,F10.3,A)') " Maximum box volume: ", max_box_volume, " cubic Angstroms"
            
            if (box_change_count > 0) then
                write(*,'(A,F10.3,A)') " Average box volume: ", &
                    total_box_volume / real(box_change_count, dp), " cubic Angstroms"
            end if
        end if
        write(*,*)
    end if
    
    call write_final_report(frame, n_atoms, end_time - start_time, atom_info, &
                          n_types, coord_numbers, atom_types)
    
    ! Report cell list statistics only if it was used
    if (use_cell_list) then
        ! Cell list statistics are reported by write_final_report
    else
        write(*,'(A)') " DIRECT CALCULATION USED"
        write(*,'(A)') " ----------------------"
        write(*,'(A)') " Cell list statistics not available when using direct calculation method"
        write(*,*)
    end if
    
    ! Cleanup analysis
    call cleanup_io()
    call cleanup_coordination()
    if (use_cell_list) then
        call cleanup_cell_list()
    end if
    call cleanup_atom_selection()
    
    ! Cleanup angle calculation if enabled
    if (calculate_angles) then
        call finalize_angle_analysis()
    end if
    
    !===============================================================================
    ! Final cleanup
    !===============================================================================
    call cleanup_config()
    
    if (allocated(coords)) deallocate(coords)
    if (allocated(atom_types)) deallocate(atom_types)
    if (allocated(elements)) deallocate(elements)
    if (allocated(include_mask)) deallocate(include_mask)
    if (allocated(atom_info)) deallocate(atom_info)
    if (allocated(type_counts)) deallocate(type_counts)
    if (allocated(explicit_pairs)) deallocate(explicit_pairs)
    
    ! Print ending timestamp
    call date_and_time(date, time, zone, values)
    write(*,'(A)') " ==========================================================================="
    write(*,'(A,A4,A1,A2,A1,A2,1X,A2,A1,A2,A1,A2)') &
              " Analysis finished at: ", date(1:4), "-", date(5:6), "-", date(7:8), &
              time(1:2), ":", time(3:4), ":", time(5:6)
    write(*,'(A)') " ==========================================================================="

contains
    !> @brief Get array of cutoff distances for all atom type pairs
    function get_cutoffs() result(cutoffs)
        real(dp), allocatable :: cutoffs(:)
        integer :: idx
        
        allocate(cutoffs(n_pairs))
        do idx = 1, n_pairs
            cutoffs(idx) = pairs(idx)%cutoff
        end do
    end function get_cutoffs
    
    !> @brief Read a trajectory frame without updating the frame counter
    !> 
    !> Helper function to read atom types without affecting main processing
    subroutine read_trajectory_frame_noupdate(coords_out, atom_types_out, elements_out, box_length_out, &
                                     frame_number, eof_out, atom_info_in, filename)
        real(dp), intent(out) :: coords_out(:,:)
        integer, intent(out) :: atom_types_out(:)
        character(len=*), intent(out) :: elements_out(:)
        real(dp), intent(out) :: box_length_out(3)
        integer, intent(inout) :: frame_number
        logical, intent(out) :: eof_out
        type(atom_type_info), intent(in) :: atom_info_in(:)
        character(len=*), intent(in) :: filename
        
        integer :: input_unit, io_stat_local, jj
        integer :: idx
        character(256) :: line
        real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
        
        ! Open file
        open(newunit=input_unit, file=filename, status='old', action='read', iostat=io_stat_local)
        if (io_stat_local /= 0) then
            eof_out = .true.
            return
        end if
        
        ! Read header - Look for timestep
        do
            read(input_unit, '(A)', iostat=io_stat_local) line
            if (io_stat_local /= 0) then
                eof_out = .true.
                close(input_unit)
                return
            end if
            
            if (index(line, 'ITEM: NUMBER OF ATOMS') > 0) exit
        end do
        
        ! Read atom count
        read(input_unit, *, iostat=io_stat_local) jj
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        ! Read box bounds if using trajectory box
        read(input_unit, '(A)', iostat=io_stat_local) line  ! ITEM: BOX BOUNDS
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        read(input_unit, *, iostat=io_stat_local) xlo, xhi
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        read(input_unit, *, iostat=io_stat_local) ylo, yhi
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        read(input_unit, *, iostat=io_stat_local) zlo, zhi
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        ! Calculate box lengths if using trajectory box dimensions
        if (use_trajectory_box) then
            box_length_out(1) = xhi - xlo
            box_length_out(2) = yhi - ylo
            box_length_out(3) = zhi - zlo
        end if
        
        ! Skip to atoms section
        read(input_unit, '(A)', iostat=io_stat_local) line  ! ITEM: ATOMS
        if (io_stat_local /= 0) then
            eof_out = .true.
            close(input_unit)
            return
        end if
        
        ! Read atom data - simplified version
        do idx = 1, min(n_atoms, size(coords_out, 1))
            read(input_unit, *, iostat=io_stat_local) jj, atom_types_out(idx), &
                coords_out(idx,1), coords_out(idx,2), coords_out(idx,3)
            if (io_stat_local /= 0) exit
            
            ! Set element name from atom_info
            if (atom_types_out(idx) >= 1 .and. atom_types_out(idx) <= size(atom_info_in)) then
                elements_out(idx) = atom_info_in(atom_types_out(idx))%name
            else
                write(elements_out(idx),'(A,I0)') 'Type', atom_types_out(idx)
            end if
        end do
        
        close(input_unit)
        eof_out = .false.
    end subroutine read_trajectory_frame_noupdate

end program coordination_analysis
