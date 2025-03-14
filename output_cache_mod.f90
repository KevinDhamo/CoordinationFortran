!> @brief Module for caching output data before writing to disk
!>
!> This module provides functionality to cache output data in memory
!> before writing it to disk, which can significantly improve I/O performance
!> by reducing the number of disk writes. It dynamically resizes the cache
!> as needed to accommodate large outputs.
module output_cache_mod
    use types_mod
    use config_mod, only: VERBOSE
    use error_mod
    implicit none
    private

    !===============================================================================
    ! Public procedures
    !===============================================================================
    public :: initialize_output_cache
    public :: cache_line
    public :: flush_cache_to_file
    public :: cleanup_output_cache
    public :: cache_header

    !===============================================================================
    ! Module variables
    !===============================================================================
    !> Path to the output file for cache flushing
    character(len=:), allocatable :: output_filename
    
    !> Flag indicating if cache has been initialized
    logical :: cache_initialized = .false.
    
    !===============================================================================
    ! Cache parameters
    !===============================================================================
    !> Initial size of the cache (number of lines)
    integer, parameter :: INITIAL_CACHE_SIZE = 1000000
    
    !> Maximum length of each cached line
    integer, parameter :: MAX_LINE_LENGTH = 2048
    
    !===============================================================================
    ! Cache storage
    !===============================================================================
    !> Array of cached output lines
    character(len=MAX_LINE_LENGTH), allocatable :: line_cache(:)
    
    !> Current number of lines in the cache
    integer :: cache_size = 0
    
    !> Current capacity of the cache (max lines without resize)
    integer :: cache_capacity = 0

contains
    !> @brief Initialize the output cache system
    !>
    !> Sets up the output cache with the specified filename and
    !> initial capacity. The cache will be flushed to this filename
    !> when flush_cache_to_file is called.
    !>
    !> @param[in] filename Path to the output file
    subroutine initialize_output_cache(filename)
        character(len=*), intent(in) :: filename
        integer :: alloc_stat
        
        ! Store the filename for later use
        if (allocated(output_filename)) deallocate(output_filename)
        allocate(character(len=len_trim(filename)) :: output_filename, stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate output filename", ERR_ALLOCATION)
            return
        end if
        output_filename = trim(filename)
        
        ! Allocate initial cache
        if (allocated(line_cache)) deallocate(line_cache)
        allocate(line_cache(INITIAL_CACHE_SIZE), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate output cache", ERR_ALLOCATION)
            return
        end if
        
        ! Initialize counters
        cache_size = 0
        cache_capacity = INITIAL_CACHE_SIZE
        cache_initialized = .true.
        
        if (VERBOSE) write(*,'(A,I0,A)') " Output caching enabled. Initial capacity: ", cache_capacity, " lines"
    end subroutine initialize_output_cache

    !> @brief Cache header information for coordination analysis output
    !>
    !> Adds standard headers to the output cache, including column format
    !> descriptions and notes about atom selection and neighbor tracking.
    !>
    !> @param[in] pairs Array of atom type pairs
    !> @param[in] n_pairs Number of atom type pairs
    !> @param[in] selection_active Flag indicating if atom selection is active
    !> @param[in] show_neighbors Flag indicating if neighbor tracking is enabled
    subroutine cache_header(pairs, n_pairs, selection_active, show_neighbors)
        type(pair_type), intent(in) :: pairs(:)
        integer, intent(in) :: n_pairs
        logical, intent(in) :: selection_active, show_neighbors
        integer :: i
        
        ! Write standard headers
        call cache_line("# Coordination Number Analysis")
        call cache_line("# Column format:")
        call cache_line("#   1. Atom ID")
        call cache_line("#   2. Atom Type")
        
        ! Write pair column headers - with no line breaks
        do i = 1, n_pairs
            write(line_cache(cache_size+1), '(A,I0,A,I0,A,I0)') "#   ", i+2, &
                ". CN(Type ", pairs(i)%type1_id, "-", pairs(i)%type2_id, ")"
            cache_size = cache_size + 1
        end do
        
        ! Add notes
        if (selection_active) then
            call cache_line("# Note: Only selected atoms are included in this output")
        end if
        
        if (show_neighbors) then
            call cache_line("# Note: Neighbor atom IDs are shown in parentheses after coordination numbers")
        end if
        
        call cache_line("# Frame data follows")
    end subroutine cache_header

    !> @brief Add a line to the output cache
    !>
    !> Adds a single line of text to the cache. Automatically resizes
    !> the cache if it's full. Each line is stored as a separate entry
    !> in the cache array.
    !>
    !> @param[in] line Text line to add to the cache
    subroutine cache_line(line)
        character(len=*), intent(in) :: line
        integer :: new_capacity, alloc_stat
        character(len=MAX_LINE_LENGTH), allocatable :: temp_cache(:)
        
        if (.not. cache_initialized) then
            call handle_error("Output cache not initialized", ERR_INVALID_PARAM, fatal=.false.)
            return
        end if
        
        ! Check if cache needs to grow
        if (cache_size >= cache_capacity) then
            ! Double the capacity
            new_capacity = cache_capacity * 2
            
            ! Allocate temporary larger cache
            allocate(temp_cache(new_capacity), stat=alloc_stat)
            if (alloc_stat /= 0) then
                call handle_error("Failed to resize output cache", ERR_ALLOCATION)
                return
            end if
            
            ! Copy existing content
            temp_cache(1:cache_size) = line_cache(1:cache_size)
            
            ! Swap old and new caches
            call move_alloc(from=temp_cache, to=line_cache)
            cache_capacity = new_capacity
            
            if (VERBOSE) write(*,'(A,I0,A)') " Output cache resized to ", cache_capacity, " lines"
        end if
        
        ! Add the new line to the cache
        cache_size = cache_size + 1
        line_cache(cache_size) = trim(line)
    end subroutine cache_line

    !> @brief Flush cached output to the output file
    !>
    !> Writes all cached lines to the output file and resets the cache.
    !> This is typically called at the end of analysis or when the cache
    !> needs to be emptied (e.g., to free memory).
    subroutine flush_cache_to_file()
        integer :: io_stat, i, output_unit
        
        if (.not. cache_initialized) then
            call handle_error("Output cache not initialized", ERR_INVALID_PARAM, fatal=.false.)
            return
        end if
        
        if (cache_size == 0) then
            if (VERBOSE) write(*,'(A)') " No data to write - output cache is empty"
            return
        end if
        
        ! Open output file
        open(newunit=output_unit, file=output_filename, status='replace', &
             action='write', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open output file: " // trim(output_filename), ERR_FILE_IO)
            return
        end if
        
        ! Write all cached lines to file
        do i = 1, cache_size
            write(output_unit, '(A)') trim(line_cache(i))
        end do
        
        ! Close file
        close(output_unit)
        
        if (VERBOSE) write(*,'(A,I0,A)') " Wrote ", cache_size, " lines to output file"
    end subroutine flush_cache_to_file

    !> @brief Clean up output cache resources
    !>
    !> Deallocates all memory used by the output cache system.
    !> This should be called at program termination to free resources.
    subroutine cleanup_output_cache()
        if (allocated(line_cache)) deallocate(line_cache)
        if (allocated(output_filename)) deallocate(output_filename)
        cache_size = 0
        cache_capacity = 0
        cache_initialized = .false.
    end subroutine cleanup_output_cache

end module output_cache_mod