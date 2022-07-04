module test_pr
    use constants
    use cubic_eos, only: pr
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_pr

contains

    subroutine collect_pr(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("peng robinson delta1", test_del1_definition), &
        new_unittest("peng robinson attractive", test_attractive_parameter) &
        ]
    end subroutine collect_pr

    subroutine test_attractive_parameter(error)
        type(error_type), allocatable, intent(out) :: error
        type(pr) :: compound1, compound2, compound3, compounds(3)
        real(wp) :: a_values(3), real_values(3)
        integer :: i

        real_values = [2282.1884, 17672.3093, 3373.5233]

        compound1 = pr("methane", 1, 2, 3, 4, 5, 6)
        compound2 = pr("nitrogen", 7, 8, 9, 10, 11, 12)
        compound3 = pr("ethane", 9, 5, 2, 15, 1, 2)

        compounds = [compound1, compound2, compound3]

        do i = 1, 3
            call compounds(i)%a_t(250._wp)
            a_values(i) = compounds(i)%a
        end do
        call check(error, maxval(abs(a_values - real_values)) > ERRMAX)
        if (allocated(error)) return
    end subroutine test_attractive_parameter

    subroutine test_del1_definition(error)
        type(error_type), allocatable, intent(out) :: error
        type(pr) :: compound

        compound = pr("methane", 1, 2, 3, 4, 5, 6)

        call check(error, 1._wp == compound%del1)
        if (allocated(error)) return

    end subroutine test_del1_definition
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

    testsuites = [&
        new_testsuite("test_pr", collect_pr) &
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
