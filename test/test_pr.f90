module test_pr
    use constants
    use cubic_eos, only: pr
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_pr

contains

    subroutine collect_pr(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("valid", test_valid) &
        ]
    end subroutine collect_pr

    subroutine test_valid(error)
        type(error_type), allocatable, intent(out) :: error
        type(pr) :: compound

        compound = pr("methane", 1, 2, 3, 4, 5, 6)

        call check(error, abs(compound%del1 - 1_wp) > ERRMAX)
        if (allocated(error)) return

    end subroutine test_valid

end module

program main
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_pr, only: collect_pr

    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("test_pr", collect_pr)]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program main
