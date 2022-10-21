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

program main
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_properties, only: collect_scalar

    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [&
        new_testsuite("collect_scalar", collect_scalar) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program main
