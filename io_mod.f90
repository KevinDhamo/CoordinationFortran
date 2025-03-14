!> @brief Module for integrating and coordinating all I/O operations
!>
!> This module acts as a facade/wrapper for the various I/O modules,
!> providing a simplified interface for initializing, using, and cleaning up
!> I/O resources. It re-exports key functions from trajectory_io_mod,
!> data_file_io_mod, output_io_mod, and progress_mod to present a unified
!> interface to the rest of the program.
module io_mod
    use types_mod
    use config_mod, only: VERBOSE, INPUT_FILE, OUTPUT_FILE  ! Add VERBOSE to the import list
    use trajectory_io_mod, only: init_trajectory_reader, read_trajectory_frame, &
                                count_trajectory_frames, cleanup_trajectory_reader, &
                                seek_to_frame, set_user_header_line
    use data_file_io_mod, only: read_data_file, read_data_file_dims_only
    use output_io_mod, only: init_output, write_output_frame, &
                           write_final_report, cleanup_output, disable_output
    use progress_mod, only: init_progress, update_progress, &
                          finalize_progress
    implicit none
    private

    !===============================================================================
    ! Re-exported procedures from other I/O modules
    !===============================================================================
    
    !> @name Data file I/O procedures
    !> @{
    
    !> Read atom type information from LAMMPS data file
    public :: read_data_file
    
    !> Read only dimensions from LAMMPS data file
    public :: read_data_file_dims_only
    !> @}
    
    !> @name Trajectory I/O procedures
    !> @{
    
    !> Read a single frame from trajectory file
    public :: read_trajectory_frame
    
    !> Seek to a specific frame in trajectory file
    public :: seek_to_frame
    
    !> Count the number of frames in trajectory file
    public :: count_trajectory_frames
    !> @}
    
    !> @name Output I/O procedures
    !> @{
    
    !> Write coordination data for a single frame
    public :: write_output_frame
    
    !> Write final analysis report
    public :: write_final_report
    
    !> Enable/disable output file writing
    public :: disable_output
    !> @}
    
    !> @name Progress tracking procedures
    !> @{
    
    !> Update progress display
    public :: update_progress
    !> @}
    
    !===============================================================================
    ! Module-specific procedures
    !===============================================================================
    
    !> Initialize all I/O systems
    public :: initialize_io
    
    !> Set trajectory header format
    public :: set_trajectory_header
    
    !> Clean up all I/O resources
    public :: cleanup_io

contains
    !> @brief Initialize all I/O systems
    !>
    !> Sets up trajectory reader, output file, and progress tracking.
    !> Optionally seeks to a specific starting frame if specified.
    !>
    !> @param[in] trajectory_file Path to input trajectory file
    !> @param[in] output_file Path to output data file
    !> @param[in] n_frames Total number of frames to process
    !> @param[in] start_at_frame Optional frame to start processing at (default: 0)
    subroutine initialize_io(trajectory_file, output_file, n_frames, start_at_frame)
        character(len=*), intent(in) :: trajectory_file, output_file
        integer, intent(in) :: n_frames
        integer, intent(in), optional :: start_at_frame
        integer :: frame_to_seek = 0
        
        ! Initialize trajectory reader
        call init_trajectory_reader(trajectory_file)
        
        ! Initialize output file
        call init_output(output_file)
        
        ! Initialize progress bar
        call init_progress(n_frames)
        
        ! Seek to starting frame if specified
        if (present(start_at_frame)) then
            frame_to_seek = start_at_frame
            if (frame_to_seek > 0) then
                call seek_to_frame(frame_to_seek)
                if (VERBOSE) then
                    write(*,'(A,I0)') ' Seeking to frame ', frame_to_seek
                end if
            end if
        end if
    end subroutine initialize_io

    !> @brief Set trajectory header format
    !>
    !> Sets a custom header format for trajectory file parsing.
    !> Used when the trajectory file has non-standard column formats.
    !>
    !> @param[in] header_format Custom header format string
    subroutine set_trajectory_header(header_format)
        character(len=*), intent(in) :: header_format
        
        if (len_trim(header_format) > 0) then
            call set_user_header_line(header_format)
        end if
    end subroutine set_trajectory_header

    !> @brief Clean up all I/O resources
    !>
    !> Properly closes all files and deallocates resources used
    !> by the various I/O modules.
    subroutine cleanup_io()
        call cleanup_trajectory_reader()
        call cleanup_output()
        call finalize_progress()
    end subroutine cleanup_io

end module io_mod