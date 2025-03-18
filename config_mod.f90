!> @brief Module for configuration and program settings management
!>
!> This module handles all configuration-related functionality,
!> including reading settings from configuration files, managing
!> default parameters, and parsing simulation setup files. It provides
!> a centralized system for configuration management across the program.
module config_mod
    use types_mod
    use error_mod
    implicit none
    private

    ! Public procedures
    public :: initialize_config
    public :: read_config
    public :: read_setup_file
    public :: get_config_value
    public :: cleanup_config

    !===============================================================================
    ! File names and paths
    !===============================================================================
    !> Path to the input trajectory file
    character(len=*), parameter, public :: INPUT_FILE = 'trajectory.dump'
    
    !> Path to the output coordination data file
    character(len=*), parameter, public :: OUTPUT_FILE = 'coordination_numbers.dat'
    
    !> Path to the input LAMMPS data file
    character(len=*), parameter, public :: DATA_FILE = 'lammps.data'
    
    !> Path to the general configuration file
    character(len=*), parameter, public :: CONFIG_FILE = 'analysis_config.txt'
    
    !> Path to the simulation setup file
    character(len=*), parameter, public :: SETUP_FILE = 'setup.txt'

    !===============================================================================
    ! Default parameters
    !===============================================================================
    !> Default cutoff distance for atom pair interactions
    real(dp), parameter, public :: DEFAULT_CUTOFF = 3.0_dp
    
    !> Width of the progress bar display
    integer, parameter, public :: PROGRESS_BAR_WIDTH = 50
    
    !> Whether to display verbose output
    logical, parameter, public :: VERBOSE = .false. ! Reduced output verbosity

    !===============================================================================
    ! Cell list configuration
    !===============================================================================
    !> Frequency of cell list updates (in frames)
    integer, public :: cell_update_freq = 5
    
    !> Factor to multiply cutoff by to determine cell size
    real(dp), parameter, public :: CELL_SIZE_FACTOR = 1.0_dp
    
    !> Initial maximum number of atoms per cell
    integer, parameter, public :: MAX_ATOMS_PER_CELL = 50
    
    !> Threshold of atom movement that triggers cell list rebuild
    real(dp), parameter, public :: REBUILD_THRESHOLD = 0.1_dp

    !===============================================================================
    ! Verlet list configuration (kept for compatibility but not used)
    !===============================================================================
    !> Whether to use Verlet lists for coordination calculations
    logical, public :: use_verlet_lists = .false.
    
    !> Skin distance around cutoff for Verlet lists
    real(dp), public :: verlet_skin_distance = 1.0_dp
    
    !> Initial capacity for Verlet list neighbor arrays
    integer, parameter, public :: INITIAL_VERLET_CAPACITY = 50
    
    !> Threshold of atom movement that triggers Verlet list rebuild
    real(dp), public :: verlet_rebuild_threshold = 0.5_dp
    
    !> Whether to show Verlet list statistics in output
    logical, public :: show_verlet_stats = .false.

    !===============================================================================
    ! Setup file parameters
    !===============================================================================
    !> Flag indicating if setup file was found and read
    logical, public :: setup_file_found = .false.
    
    !> Number of CPU cores requested for parallel execution
    integer, public :: requested_cores = 0
    
    !> Start processing from this frame (0 = beginning)
    integer, public :: start_frame = 0
    
    !> End processing at this frame (0 = all frames)
    integer, public :: end_frame = 0
    
    !> Original string for end frame from setup file
    character(len=32), public :: end_frame_str = "unlimited"
    
    !> Enable selective cell list rebuilding for performance
    logical, public :: selective_rebuild = .false.
    
    !> Path to the input data file specified in setup
    character(len=256), public :: input_data_file = ''
    
    !> Path to the input trajectory file specified in setup
    character(len=256), public :: input_trajectory = ''
    
    !> Path to the output data file specified in setup
    character(len=256), public :: output_data_file = ''
    
    !===============================================================================
    ! Enhanced functionality parameters
    !===============================================================================
    !> User-provided header format for trajectory parsing
    character(len=512), public :: trajectory_header_format = ''
    
    !> Atom selection string (e.g., "1-100,150,200-300")
    character(len=1024), public :: selected_atoms = ''
    
    !> Flag to show atom neighbor IDs in output
    logical, public :: show_atom_neighbors = .false.

    !> Flag to use box dimensions from trajectory file instead of data file
    logical, public :: use_trajectory_box = .false.

    !===============================================================================
    ! Type definition for configuration key-value pairs
    !===============================================================================
    !> Type for storing configuration key-value pairs
    type :: config_value_type
        character(len=32) :: key = ''
        character(len=128) :: value = ''
    end type config_value_type

    ! Module variables
    type(config_value_type), allocatable :: config_data(:)

contains

    !> @brief Initialize configuration with default values
    !>
    !> Sets up the configuration system and initializes all
    !> parameters to their default values.
    subroutine initialize_config()
        integer :: n_config_items = 20  ! Initial size
        integer :: alloc_stat
        
        ! Allocate configuration array
        if (allocated(config_data)) deallocate(config_data)
        allocate(config_data(n_config_items), stat=alloc_stat)
        call check_allocation(alloc_stat, "config_data")
        
        ! Initialize all entries to empty strings
        config_data%key = ''
        config_data%value = ''
        
        ! Set defaults for setup parameters
        setup_file_found = .false.
        requested_cores = 0
        start_frame = 0
        end_frame = 0
        end_frame_str = "unlimited"
        cell_update_freq = 5  ! Default value of 5, can be overridden in setup.txt
        selective_rebuild = .false.
        input_data_file = DATA_FILE
        input_trajectory = INPUT_FILE
        output_data_file = OUTPUT_FILE
        trajectory_header_format = ''
        selected_atoms = ''
        show_atom_neighbors = .false.
        use_trajectory_box = .false.
        
        ! Initialize Verlet list parameters (disabled by default)
        use_verlet_lists = .false.
        verlet_skin_distance = 1.0_dp
        verlet_rebuild_threshold = 0.5_dp
        show_verlet_stats = .false.
    end subroutine initialize_config

    !> @brief Read configuration from the configuration file
    !>
    !> Reads key-value pairs from the configuration file and
    !> stores them in the configuration system.
    subroutine read_config()
        integer :: io_stat, unit
        character(len=256) :: line
        character(len=32) :: key
        character(len=128) :: value
        logical :: file_exists
        
        ! Check if configuration file exists
        inquire(file=CONFIG_FILE, exist=file_exists)
        if (.not. file_exists) return
        
        ! Open configuration file
        open(newunit=unit, file=CONFIG_FILE, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open configuration file", ERR_FILE_IO)
            return
        end if
        
        ! Read and parse each line of the file
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse key-value pair
            call parse_config_line(line, key, value)
            if (len_trim(key) > 0) then
                call store_config_value(key, value)
            end if
        end do
        
        close(unit)
    end subroutine read_config
    
    !> @brief Parse cutoff key to extract type indices
    !>
    !> Parses a pair cutoff key string (e.g., "PAIR_CUTOFF_1_2")
    !> to extract the atom type indices.
    !>
    !> @param[in] key_str Key string to parse
    !> @param[out] type1 First atom type index
    !> @param[out] type2 Second atom type index
    !> @param[out] success Flag indicating if parsing was successful
    subroutine parse_cutoff_key(key_str, type1, type2, success)
        character(len=*), intent(in) :: key_str
        integer, intent(out) :: type1, type2
        logical, intent(out) :: success
        
        integer :: underscore_pos, io_stat
        character(len=32) :: prefix, type1_str, type2_str
        
        ! Initialize outputs
        success = .false.
        type1 = -1
        type2 = -1
        
        ! Check if it's a cutoff key
        prefix = "PAIR_CUTOFF_"
        if (len_trim(key_str) <= len_trim(prefix)) return
        if (key_str(1:len_trim(prefix)) /= prefix) return
        
        ! Extract the part after the prefix
        type1_str = key_str(len_trim(prefix)+1:)
        
        ! Find the underscore
        underscore_pos = index(type1_str, "_")
        if (underscore_pos <= 1) return
        
        ! Split at the underscore
        type2_str = type1_str(underscore_pos+1:)
        type1_str = type1_str(1:underscore_pos-1)
        
        ! Convert to integers
        read(type1_str, *, iostat=io_stat) type1
        if (io_stat /= 0) return
        
        read(type2_str, *, iostat=io_stat) type2
        if (io_stat /= 0) return
        
        success = .true.
    end subroutine parse_cutoff_key
 
    !> @brief Read simulation setup from the setup file
    !>
    !> Reads simulation parameters, atom types, and pair cutoffs
    !> from the setup file. Creates atom_info and pairs arrays.
    !>
    !> @param[out] atom_info Array of atom type information
    !> @param[out] n_types Number of atom types
    !> @param[out] pairs Array of pair interaction parameters
    !> @param[out] n_pairs Number of atom type pairs
    !> @return .true. if setup file was read successfully, .false. otherwise
    function read_setup_file(atom_info, n_types, pairs, n_pairs) result(success)
        type(atom_type_info), allocatable, intent(out) :: atom_info(:)
        integer, intent(out) :: n_types
        type(pair_type), allocatable, intent(out) :: pairs(:)
        integer, intent(out) :: n_pairs
        logical :: success
        
        integer :: io_stat, unit, i, j, pair_idx
        character(len=512) :: line, value_str
        character(len=32) :: key
        logical :: file_exists
        logical :: atom_types_found = .false.
        character(len=8), allocatable :: type_names(:)
        logical, allocatable :: type_include(:)
        real(dp), allocatable :: type_masses(:)
        character(len=256) :: comment
        integer :: comment_pos
        logical :: success_parse
        
        ! Initialize return value
        success = .false.
        
        ! Check if setup file exists
        inquire(file=SETUP_FILE, exist=file_exists)
        if (.not. file_exists) then
            if (VERBOSE) write(*,*) "Setup file not found, using interactive mode"
            return
        end if
        
        ! Open setup file
        open(newunit=unit, file=SETUP_FILE, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open setup file", ERR_FILE_IO, fatal=.false.)
            return
        end if
        
        ! First pass: count atom types and find basic parameters
        n_types = 0
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            key = adjustl(line)
            value_str = ''
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for comments
                comment_pos = index(value_str, '#')
                if (comment_pos > 0) then
                    comment = adjustl(value_str(comment_pos+1:))
                    value_str = adjustl(value_str(:comment_pos-1))
                end if
            end if
            
            ! Process key-value pairs
            select case(trim(key))
                case('ATOM_TYPES')
                    read(value_str, *, iostat=io_stat) n_types
                    if (io_stat == 0 .and. n_types > 0) then
                        atom_types_found = .true.
                        allocate(type_names(n_types), type_include(n_types), type_masses(n_types))
                        type_names = ''
                        type_include = .true.
                        type_masses = 0.0_dp
                    end if
                
                case('Cores')
                    read(value_str, *, iostat=io_stat) requested_cores
                
                case('START_FRAME')
                    read(value_str, *, iostat=io_stat) start_frame
                    if (io_stat /= 0) start_frame = 0
                
                case('END_FRAME')
                    if (trim(value_str) == 'unlimited' .or. &
                        trim(value_str) == 'all') then
                        end_frame = 0  ! Indicating unlimited
                        end_frame_str = trim(value_str)  ! Store original string
                    else
                        read(value_str, *, iostat=io_stat) end_frame
                        if (io_stat /= 0) end_frame = 0
                        ! Store numeric value as string
                        if (io_stat == 0) then
                            write(end_frame_str, '(I0)') end_frame
                        else
                            end_frame_str = "unlimited"
                        end if
                    end if
                    
                case('CELL_UPDATE_FREQ')
                    read(value_str, *, iostat=io_stat) cell_update_freq
                    if (io_stat /= 0) cell_update_freq = 5  ! Default to 5 if invalid input
                
                case('SELECTIVE_REBUILD')
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        selective_rebuild = .true.
                    else
                        selective_rebuild = .false.
                    end if
                
                case('INPUT_DATA')
                    input_data_file = trim(value_str)
                
                case('INPUT_TRAJECTORY')
                    input_trajectory = trim(value_str)
                
                case('OUTPUT_FILE')
                    output_data_file = trim(value_str)
                    
                ! New parameters
                case('TRAJECTORY_HEADER')
                    trajectory_header_format = trim(value_str)
                
                case('SELECTED_ATOMS')
                    selected_atoms = trim(value_str)
                    
                case('SHOW_NEIGHBORS')
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        show_atom_neighbors = .true.
                    else
                        show_atom_neighbors = .false.
                    end if
                
                case('USE_TRAJECTORY_BOX')
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        use_trajectory_box = .true.
                    else
                        use_trajectory_box = .false.
                    end if
                    
                ! Verlet list parameters (kept for compatibility but ignored)
                case('USE_VERLET_LISTS')
                    ! Keep parsing but ignore the value - we always use false
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        ! Silently ignore - we're not using Verlet lists
                    else
                        ! This is the default anyway
                    end if
                    
                case('VERLET_SKIN_DISTANCE')
                    read(value_str, *, iostat=io_stat) verlet_skin_distance
                    if (io_stat /= 0) verlet_skin_distance = 1.0_dp
                    if (verlet_skin_distance <= 0.0_dp) verlet_skin_distance = 1.0_dp
                    
                case('VERLET_REBUILD_THRESHOLD')
                    read(value_str, *, iostat=io_stat) verlet_rebuild_threshold
                    if (io_stat /= 0) verlet_rebuild_threshold = 0.5_dp
                    if (verlet_rebuild_threshold <= 0.0_dp .or. &
                        verlet_rebuild_threshold >= 1.0_dp) then
                        verlet_rebuild_threshold = 0.5_dp
                    end if
                    
                case('SHOW_VERLET_STATS')
                    if (trim(value_str) == 'yes' .or. &
                        trim(value_str) == 'true' .or. &
                        trim(value_str) == '1') then
                        show_verlet_stats = .true.
                    else
                        show_verlet_stats = .false.
                    end if
            end select
        end do
        
        ! Verify that we found atom types
        if (.not. atom_types_found .or. n_types <= 0) then
            if (VERBOSE) write(*,*) "No valid ATOM_TYPES found in setup file"
            close(unit)
            return
        end if
        
        ! Second pass: read atom types and pair cutoffs
        rewind(unit)
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            key = adjustl(line)
            value_str = ''
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for ignore comment
                comment_pos = index(value_str, '#')
                if (comment_pos > 0) then
                    comment = adjustl(value_str(comment_pos+1:))
                    value_str = adjustl(value_str(:comment_pos-1))
                    
                    ! Check if it's an ignore directive
                    if (index(comment, 'ignore') > 0) then
                        ! Extract the type number from TYPE_X
                        if (key(1:5) == 'TYPE_') then
                            read(key(6:), *, iostat=io_stat) i
                            if (io_stat == 0 .and. i > 0 .and. i <= n_types) then
                                type_include(i) = .false.
                                if (VERBOSE) write(*,*) "Ignoring atom type:", i
                            end if
                        end if
                    end if
                end if
            end if
            
            ! Process atom types
            if (key(1:5) == 'TYPE_') then
                read(key(6:), *, iostat=io_stat) i
                if (io_stat == 0 .and. i > 0 .and. i <= n_types) then
                    type_names(i) = trim(adjustl(value_str))
                end if
            end if
        end do
        
        ! Calculate number of unique pairs
        n_pairs = 0
        do i = 1, n_types
            if (.not. type_include(i)) cycle
            do j = i, n_types
                if (.not. type_include(j)) cycle
                n_pairs = n_pairs + 1
            end do
        end do
        
        ! Allocate array for atom types
        if (allocated(atom_info)) deallocate(atom_info)
        allocate(atom_info(n_types), stat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to allocate atom_info array", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize atom info
        do i = 1, n_types
            atom_info(i)%type_id = i
            atom_info(i)%name = type_names(i)
            atom_info(i)%mass = type_masses(i)
            atom_info(i)%include = type_include(i)
        end do
        
        ! Allocate pairs array
        if (allocated(pairs)) deallocate(pairs)
        allocate(pairs(n_pairs), stat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to allocate pairs array", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize pair indices first
        pair_idx = 1
        do i = 1, n_types
            if (.not. type_include(i)) cycle
            do j = i, n_types
                if (.not. type_include(j)) cycle
                
                ! Set default cutoff
                pairs(pair_idx)%type1_id = i
                pairs(pair_idx)%type2_id = j
                pairs(pair_idx)%cutoff = DEFAULT_CUTOFF
                pairs(pair_idx)%cutoff_sq = DEFAULT_CUTOFF**2
                
                pair_idx = pair_idx + 1
            end do
        end do
        
        ! Read cutoff values
        rewind(unit)
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse line
            if (index(line, '=') > 0) then
                key = adjustl(line(:index(line, '=')-1))
                value_str = adjustl(line(index(line, '=')+1:))
                
                ! Check for PAIR_CUTOFF entries
                if (key(1:12) == 'PAIR_CUTOFF_') then
                    ! Parse the cutoff key to extract type indices
                    call parse_cutoff_key(key, i, j, success_parse)
                    
                    if (success_parse .and. i > 0 .and. i <= n_types .and. j > 0 .and. j <= n_types) then
                        ! Find the pair index
                        do pair_idx = 1, n_pairs
                            ! Ensure i <= j for consistent lookup
                            if ((pairs(pair_idx)%type1_id == min(i,j)) .and. &
                                (pairs(pair_idx)%type2_id == max(i,j))) then
                                
                                read(value_str, *, iostat=io_stat) pairs(pair_idx)%cutoff
                                if (io_stat == 0) then
                                    pairs(pair_idx)%cutoff_sq = pairs(pair_idx)%cutoff**2
                                end if
                                exit
                            end if
                        end do
                    end if
                end if
            end if
        end do
        
        close(unit)
        
        ! Cleanup temporary arrays
        if (allocated(type_names)) deallocate(type_names)
        if (allocated(type_include)) deallocate(type_include)
        if (allocated(type_masses)) deallocate(type_masses)
        
        setup_file_found = .true.
        success = .true.
    end function read_setup_file

    !> @brief Parse a configuration line into key and value
    !>
    !> Parses a line of the form "key = value" into separate key and value strings.
    !>
    !> @param[in] line Line to parse
    !> @param[out] key Extracted key
    !> @param[out] value Extracted value
    subroutine parse_config_line(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key
        character(len=*), intent(out) :: value
        integer :: eq_pos
        
        ! Find equals sign in the line
        eq_pos = index(line, '=')
        if (eq_pos > 0) then
            key = adjustl(line(:eq_pos-1))
            value = adjustl(line(eq_pos+1:))
        else
            key = ''
            value = ''
        end if
    end subroutine parse_config_line

    !> @brief Store a configuration key-value pair
    !>
    !> Stores a configuration key-value pair in the config_data array.
    !> If the key already exists, its value is updated.
    !>
    !> @param[in] key Configuration key
    !> @param[in] value Configuration value
    subroutine store_config_value(key, value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value
        integer :: i, empty_slot
        
        ! Find existing key or first empty slot
        empty_slot = -1
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                config_data(i)%value = value
                return
            else if (empty_slot == -1 .and. len_trim(config_data(i)%key) == 0) then
                empty_slot = i
            end if
        end do
        
        ! Store in empty slot if found
        if (empty_slot > 0) then
            config_data(empty_slot)%key = key
            config_data(empty_slot)%value = value
        else
            call handle_error("Configuration storage full", ERR_INVALID_PARAM, fatal=.false.)
        end if
    end subroutine store_config_value

    !> @brief Get a configuration value
    !>
    !> Retrieves a configuration value for the specified key.
    !> If the key is not found, returns the provided default value.
    !>
    !> @param[in] key Configuration key to look up
    !> @param[in] default Default value to return if key not found
    !> @return Configuration value or default
    function get_config_value(key, default) result(value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: default
        character(len=128) :: value
        integer :: i
        
        ! Set default value initially
        value = default
        
        ! Return default if config_data not allocated
        if (.not. allocated(config_data)) return
        
        ! Search for key in config_data
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                value = config_data(i)%value
                return
            end if
        end do
    end function get_config_value

    !> @brief Clean up configuration resources
    !>
    !> Deallocates the configuration data array.
    subroutine cleanup_config()
        if (allocated(config_data)) deallocate(config_data)
    end subroutine cleanup_config

end module config_mod
