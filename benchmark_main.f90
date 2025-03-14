!> @brief Standalone program for benchmarking coordination analysis performance
!>
!> This is a separate main program that implements benchmarking functionality.
!> It is essentially a modified version of main.f90 with built-in benchmarking logic.
program benchmark_main
    use types_mod
    use config_mod
    use cell_list_mod
    use coordination_mod
    use io_mod
    use error_mod
    use atom_selection_mod
    use omp_lib
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
    
    ! Benchmark variables
    integer :: max_threads                   ! Maximum available OpenMP threads
    integer :: bench_frame                   ! Frame counter for current benchmark
    real(dp) :: bench_start_time             ! Start time for current benchmark
    real(dp) :: bench_end_time               ! End time for current benchmark
    real(dp) :: bench_total_time             ! Total time for benchmark run
    real(dp) :: frames_per_second            ! Performance metric (frames/second)
    real(dp) :: single_thread_fps = 0.0_dp   ! Performance of single-thread run (for speedup)
    integer :: bench_unit                    ! File unit for benchmark results
    character(len=*), parameter :: BENCHMARK_OUTPUT = 'benchmark_results.txt'
    logical :: benchmark_file_open = .false.
    real(dp) :: speedup, efficiency
    
    ! Store original pairs for benchmark runs
    type(pair_type), allocatable :: benchmark_pairs(:)
    integer :: benchmark_n_pairs = 0
    
    ! Benchmark time limit (seconds)
    integer, parameter :: BENCHMARK_TIME_LIMIT = 30
    
    ! Variable for type counts and percentages
    integer :: type_count, i, j, k, type1, type2, types_found, temp_id
    real(dp) :: percentage, box_volume
    integer, allocatable :: type_counts(:)
    
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
    
    !===============================================================================
    ! Read setup file
    !===============================================================================
    setup_file_read = read_setup_file(atom_info, n_types, pairs, n_pairs)
    
    !===============================================================================
    ! Make a backup copy of the pairs array from setup.txt
    !===============================================================================
    benchmark_n_pairs = n_pairs
    allocate(benchmark_pairs(n_pairs))
    do i = 1, n_pairs
        benchmark_pairs(i) = pairs(i)
    end do
        
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
    
    ! Allocate and initialize type counts array
    allocate(type_counts(n_types))
    type_counts = 0
    
    ! Read first frame to get atom types
    if (setup_file_read) then
        call read_trajectory_frame_noupdate(coords, atom_types, elements, box_length, &
                                      frame, eof, atom_info, input_trajectory)
    else
        call read_trajectory_frame_noupdate(coords, atom_types, elements, box_length, &
                                      frame, eof, atom_info, INPUT_FILE)
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
        if (atom_info(i)%include) then
            types_found = types_found + 1
            type_count = type_counts(i)
            percentage = real(type_count) / real(n_atoms) * 100.0
            
            ! Format to match ideal output
            if (i == 1) then
                write(*,'(A,I0,A,A,A,I0,A,F6.2,A)') "   Type ", i, " (", &
                    trim(atom_info(i)%name), "):             ", type_count, " atoms (", percentage, "%)"
            else if (i == 2) then
                write(*,'(A,I0,A,A,A,I0,A,F6.2,A)') "   Type ", i, " (", &
                    trim(atom_info(i)%name), "):             ", type_count, " atoms   (", percentage, "%)"
            else
                write(*,'(A,I0,A,A,A,I0,A,F6.2,A)') "   Type ", i, " (", &
                    trim(atom_info(i)%name), "):              ", type_count, " atoms  (", percentage, "%)"
            end if
        end if
    end do
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
    
    write(*,'(A,I0)') " Cell update frequency:     ", cell_update_freq
    if (selective_rebuild) then
        write(*,'(A)') " Selective cell rebuilding: True"
    else
        write(*,'(A)') " Selective cell rebuilding: False"
    end if
    
    ! Get maximum available threads for benchmarking
    max_threads = omp_get_max_threads()
    write(*,'(A,I0)') " Maximum available threads: ", max_threads
    
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
    ! Set up benchmark output file
    !===============================================================================
    write(*,'(A)') ""
    write(*,'(A)') " BENCHMARK MODE"
    write(*,'(A)') " --------------"
    write(*,'(A)') " Thread scaling performance test"
    write(*,'(A,I0,A)') " Each configuration will run for up to ", BENCHMARK_TIME_LIMIT, " seconds"
    write(*,'(A)') " Using configuration from setup.txt"
    write(*,'(A)') ""
    
    ! Open benchmark output file
    open(newunit=bench_unit, file=BENCHMARK_OUTPUT, status='replace', iostat=io_stat)
    if (io_stat /= 0) then
        call handle_error("Failed to open benchmark file", ERR_FILE_IO)
    end if
    benchmark_file_open = .true.
    
    ! Write header to benchmark file
    write(bench_unit,'(A)') "# Coordination Analysis Benchmark Results"
    write(bench_unit,'(A)') "# Using configuration from setup.txt"
    write(bench_unit,'(A)') "# Threads, Time(s), Frames/s, Frames, Speedup, Efficiency"
    
    !===========================================================================
    ! Display pair cutoffs from setup.txt
    !===========================================================================
    write(*,'(A)') " PAIR CUTOFFS FROM SETUP.TXT"
    write(*,'(A)') " --------------------------"
    do i = 1, benchmark_n_pairs
        write(*,'(A,I0,A,A,A,I0,A,A,A,F5.2,A)') &
            " Type ", benchmark_pairs(i)%type1_id, " (", &
            trim(atom_info(benchmark_pairs(i)%type1_id)%name), ") - Type ", &
            benchmark_pairs(i)%type2_id, " (", &
            trim(atom_info(benchmark_pairs(i)%type2_id)%name), "): ", &
            benchmark_pairs(i)%cutoff, " Angstrom"
    end do
    write(*,*)
    
    !===============================================================================
    ! Run benchmark with varying thread counts
    !===============================================================================
    write(*,'(A)') "=== Starting Benchmark ==="
    
    ! Disable output file writing during benchmarks
    disable_output = .true.
    
    ! Test different thread counts (powers of 2)
    num_threads = 1
    do while (num_threads <= max_threads)
        write(*,'(A,I0,A)') "Testing with ", num_threads, " threads..."
        
        ! Set OpenMP parameters - override what might be in setup.txt
        call omp_set_num_threads(num_threads)
        call omp_set_schedule(omp_sched_guided, 0)  ! Use guided scheduling
        
        ! Initialize for this run
        bench_frame = 0
        eof = .false.
        
        ! Clean up from previous run
        call cleanup_coordination()
        call cleanup_cell_list()
        call cleanup_atom_selection()
        call cleanup_io()
        
        ! Restore pairs array from our backup - THIS IS CRITICAL
        if (allocated(pairs)) deallocate(pairs)
        allocate(pairs(benchmark_n_pairs))
        do i = 1, benchmark_n_pairs
            pairs(i) = benchmark_pairs(i)
        end do
        n_pairs = benchmark_n_pairs
            
        ! Set up coordination - Initialize coordination and spatial partitioning
        call initialize_coordination(n_atoms, n_types, atom_info)
        call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
        
        ! Set up atom selection if specified
        if (len_trim(selected_atoms) > 0) then
            call initialize_atom_selection(n_atoms)
            call parse_atom_selection(selected_atoms, n_atoms)
        end if
        
        ! Initialize I/O
        if (setup_file_read) then
            call initialize_io(input_trajectory, output_data_file, total_frames, start_frame)
            bench_frame = start_frame
        else
            call initialize_io(INPUT_FILE, OUTPUT_FILE, total_frames)
        end if
        
        ! Start timing
        bench_start_time = omp_get_wtime()
            
        ! Process trajectory frames for the current thread count
        do while (.not. eof)
            call read_trajectory_frame(coords, atom_types, elements, box_length, &
                                     bench_frame, eof, atom_info)
            if (eof) exit
                
            call calculate_coordination(coords, atom_types, n_atoms, box_length, &
                                      bench_frame, include_mask)
            
            ! Update progress (with reduced frequency for benchmark)
            if (mod(bench_frame, 10) == 0) then
                write(*,'(a1,1x,A,I0,A,I0,A)',advance='no') achar(13), &
                    'Processing frame ', bench_frame, ' of ', total_frames, '    '
            end if
            
            ! Check if we've reached time limit and processed enough frames
            if ((omp_get_wtime() - bench_start_time) >= BENCHMARK_TIME_LIMIT .and. bench_frame >= 10) exit
            
            ! If we have frame selection in setup.txt, respect it
            if (end_frame > 0 .and. bench_frame >= end_frame) exit
        end do
            
        ! Clear progress line
        write(*,'(a1,A)',advance='yes') achar(13), repeat(' ', 80)
            
        ! Calculate results for this thread count
        bench_end_time = omp_get_wtime()
        bench_total_time = bench_end_time - bench_start_time
            
        ! Calculate frames per second
        if (bench_total_time > 0.0) then
            frames_per_second = real(bench_frame) / bench_total_time
        else
            frames_per_second = 0.0
        end if
            
        ! Calculate speedup and efficiency
        if (num_threads == 1) then
            single_thread_fps = frames_per_second
            speedup = 1.0_dp
            efficiency = 1.0_dp
        else if (single_thread_fps > 0.0_dp) then
            speedup = frames_per_second / single_thread_fps
            efficiency = speedup / real(num_threads, dp)
        else
            speedup = 0.0_dp
            efficiency = 0.0_dp
        end if
            
        ! Write results to screen
        write(*,'(A,I0,A)') "Results for ", num_threads, " threads:"
        write(*,'(A,F10.3,A)') "  Time: ", bench_total_time, " seconds"
        write(*,'(A,I0)') "  Frames processed: ", bench_frame
        write(*,'(A,F10.3,A)') "  Performance: ", frames_per_second, " frames/second"
        write(*,'(A,F10.3)') "  Speedup: ", speedup
        write(*,'(A,F10.3,A)') "  Efficiency: ", efficiency * 100.0, "%"
        write(*,*)
            
        ! Write to benchmark file
        write(bench_unit,'(I6,5(F12.3))') num_threads, bench_total_time, &
            frames_per_second, real(bench_frame), speedup, efficiency
        call flush(bench_unit)
            
        ! Double the thread count for next iteration
        if (num_threads == max_threads) exit
        num_threads = min(num_threads * 2, max_threads)
    end do
    
    !===============================================================================
    ! Final cleanup
    !===============================================================================
    write(*,'(A)') "Benchmark complete. Results written to: " // BENCHMARK_OUTPUT
    
    if (benchmark_file_open) then
        close(bench_unit)
        benchmark_file_open = .false.
    end if
    
    if (allocated(benchmark_pairs)) deallocate(benchmark_pairs)
    
    call cleanup_io()
    call cleanup_coordination()
    call cleanup_cell_list()
    call cleanup_atom_selection()
    call cleanup_config()
    
    if (allocated(coords)) deallocate(coords)
    if (allocated(atom_types)) deallocate(atom_types)
    if (allocated(elements)) deallocate(elements)
    if (allocated(include_mask)) deallocate(include_mask)
    if (allocated(atom_info)) deallocate(atom_info)
    if (allocated(type_counts)) deallocate(type_counts)
    
    ! Print ending timestamp
    call date_and_time(date, time, zone, values)
    write(*,'(A)') " ==========================================================================="
    write(*,'(A,A4,A1,A2,A1,A2,1X,A2,A1,A2,A1,A2)') &
              " Benchmark finished at: ", date(1:4), "-", date(5:6), "-", date(7:8), &
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
        
        ! Open file
        open(newunit=input_unit, file=filename, status='old', action='read', iostat=io_stat_local)
        if (io_stat_local /= 0) then
            eof_out = .true.
            return
        end if
        
        ! Read header - simplified for atom type detection
        do idx = 1, 9
            read(input_unit, '(A)', iostat=io_stat_local)  ! Skip header lines
            if (io_stat_local /= 0) then
                eof_out = .true.
                return
            end if
        end do
        
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

end program benchmark_main