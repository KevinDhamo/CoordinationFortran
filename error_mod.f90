!> @brief Module for error handling and parameter validation
!>
!> This module provides centralized error handling and parameter validation
!> functions for consistent error reporting throughout the application.
!> It includes facilities for reporting errors, checking allocation status,
!> and validating numerical parameters against allowed ranges.
module error_mod
    use iso_fortran_env, only: error_unit
    implicit none
    private

    ! Public procedures
    public :: handle_error
    public :: check_allocation
    public :: validate_parameters

    !> @name Error codes
    !> @{
    
    !> Error code for memory allocation failures
    integer, parameter, public :: ERR_ALLOCATION = 1
    
    !> Error code for file I/O operations
    integer, parameter, public :: ERR_FILE_IO = 2
    
    !> Error code for invalid parameter values
    integer, parameter, public :: ERR_INVALID_PARAM = 3
    
    !> Error code for invalid type specifications
    integer, parameter, public :: ERR_INVALID_TYPE = 4
    
    !> Error code for inconsistent data structures
    integer, parameter, public :: ERR_INCONSISTENT_DATA = 5
    !> @}

contains

    !> @brief Handle an error with customizable severity
    !>
    !> Reports an error message with error code to the error output unit.
    !> Can optionally terminate the program for fatal errors.
    !>
    !> @param[in] msg Error message to display
    !> @param[in] error_code Numeric error code (see module parameters)
    !> @param[in] fatal Optional flag to indicate if error is fatal (default: .true.)
    subroutine handle_error(msg, error_code, fatal)
        character(len=*), intent(in) :: msg
        integer, intent(in) :: error_code
        logical, intent(in), optional :: fatal
        logical :: is_fatal
        
        ! Determine if error is fatal (default is true)
        is_fatal = .true.
        if (present(fatal)) is_fatal = fatal

        ! Write error message to error output unit
        write(error_unit,'(A,I0,2A)') "Error (", error_code, "): ", trim(msg)
        
        ! Terminate program if error is fatal
        if (is_fatal) then
            error stop
        end if
    end subroutine handle_error

    !> @brief Check allocation status and report errors
    !>
    !> Checks the status code returned by an allocation operation
    !> and reports an error if allocation failed.
    !>
    !> @param[in] alloc_stat Status code from allocation
    !> @param[in] array_name Name of the array for error reporting
    subroutine check_allocation(alloc_stat, array_name)
        integer, intent(in) :: alloc_stat
        character(len=*), intent(in) :: array_name
        
        if (alloc_stat /= 0) then
            call handle_error("Failed to allocate memory for " // trim(array_name), &
                            ERR_ALLOCATION)
        end if
    end subroutine check_allocation

    !> @brief Validate numerical parameters against allowed ranges
    !>
    !> Checks if a parameter value is within the specified minimum and
    !> maximum values, reporting an error if validation fails.
    !>
    !> @param[in] param_name Name of the parameter for error reporting
    !> @param[in] param_value Value of the parameter to validate
    !> @param[in] min_value Minimum allowed value
    !> @param[in] max_value Maximum allowed value
    subroutine validate_parameters(param_name, param_value, min_value, max_value)
        character(len=*), intent(in) :: param_name
        real, intent(in) :: param_value
        real, intent(in) :: min_value
        real, intent(in) :: max_value
        
        if (param_value < min_value .or. param_value > max_value) then
            write(error_unit,'(5A,2(A,G0.4))') "Parameter '", trim(param_name), &
                "' value ", "out of range [", trim(param_name), &
                "] = ", param_value, " not in [", min_value, ",", max_value, "]"
            call handle_error("Parameter validation failed", ERR_INVALID_PARAM)
        end if
    end subroutine validate_parameters

end module error_mod
