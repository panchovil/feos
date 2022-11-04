module test_properties
    use constants, only: wp, errmax
    use properties
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    !private
    public :: collect_scalar

contains

    ! =============================================================================
    !  Tests that will be run
    ! -----------------------------------------------------------------------------
    subroutine collect_scalar(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("real addition", test_real_scalar_addition) &
        ]
    end subroutine collect_scalar
    ! =============================================================================

    subroutine test_real_scalar_addition(error)
        type(error_type), allocatable, intent(out) :: error
        type(scalar_property) :: scalar

        scalar%val = 2.0d0
        call check(error, abs(2.d0 + scalar) - 4.d0 <= errmax)
        call check(error, abs(scalar + 2.d0) - 4.d0 <= errmax)

        if (allocated(error)) return
    end subroutine test_real_scalar_addition

end module test_properties
