!module test_properties
!    use constants, only: wp, errmax
!    use properties ! , only: scalar_property
!    use testdrive, only : new_unittest, unittest_type, error_type, check
!
!    implicit none
!
!    !private
!    public :: collect_scalar
!
!contains
!
!    ! =============================================================================
!    !  Tests that will be run
!    ! -----------------------------------------------------------------------------
!    subroutine collect_scalar(testsuite)
!      type(unittest_type), allocatable, intent(out) :: testsuite(:)
!
!      testsuite = [ &
!        new_unittest("real addition", test_real_scalar_addition) &
!        ]
!    end subroutine collect_scalar
!    ! =============================================================================
!
!    subroutine test_real_scalar_addition(error)
!        type(error_type), allocatable, intent(out) :: error
!        type(scalar_property(2)) :: scalar
!
!        scalar = get_scalar_type()
!
!        call check(error, (2.d0 + scalar) - 4.d0 <= errmax)
!
!        if (allocated(error)) return
!    end subroutine test_real_scalar_addition
!
!    function get_scalar_type() result(scalar_type)
!        type(scalar_property(2)) :: scalar_type
!        
!        scalar_type = scalar_property(&
!            n=2, &
!            val=2, &
!            dt=3, &
!            dt2=4, &
!            dv=5, &
!            dv2=6, &
!            dtv=7, &
!            dn=[1, 2], &
!            dn2=reshape([1,2, 3, 4], [2,2]), &
!            dvn=[1,2], &
!            dtn=[1,2] &
!            )
!    end function get_scalar_type
!end module test_properties
!
!program main
!    use, intrinsic :: iso_fortran_env, only : error_unit
!    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
!    use test_properties, only: collect_scalar
!
!    implicit none
!
!    integer :: stat, is
!    type(testsuite_type), allocatable :: testsuites(:)
!    character(len=*), parameter :: fmt = '("#", *(1x, a))'
!    
!    stat = 0
!
!    testsuites = [&
!        new_testsuite("collect_scalar", collect_scalar) &
!        ]
!
!    do is = 1, size(testsuites)
!        write(error_unit, fmt) "Testing:", testsuites(is)%name
!        call run_testsuite(testsuites(is)%collect, error_unit, stat)
!    end do
!
!    if (stat > 0) then
!        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
!        error stop
!    end if
!end program main
