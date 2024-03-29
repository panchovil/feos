module test_pr
    use constants
    use properties
    use cubic_eos, only: PengRobinson, PR
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_pr

contains
    ! =============================================================================
    !  Tests that will be run
    ! -----------------------------------------------------------------------------
    subroutine collect_pr(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("PR: delta1", test_del1_definition), &
                    new_unittest("PR: get parameters", test_get_params), &
                    new_unittest("PR: get critical", test_get_critical), &
                    new_unittest("PR: attractive", test_attractive_parameter) &
                    ]
    end subroutine collect_pr
    ! =============================================================================

    subroutine test_get_params(error)
        ! Obtain EoS parameters from critical constants
        type(error_type), allocatable, intent(out) :: error
        type(PengRobinson) :: compound
        character(len=:), allocatable :: name
        real(wp) :: real_values(3), calc_values(3)

        name = "methane"

        compound = PR(name, tc=191.15_wp, pc=46.41_wp, w=0.0115_wp)

        call compound%get_params()
        real_values = [2.488, 0.0266, 0.3923]
        calc_values = [compound%ac, compound%b, compound%k]

        call check(error, maxval(abs(real_values - calc_values)) > errmax)

        if (allocated(error)) return
    end subroutine test_get_params

    subroutine test_get_critical(error)
        type(error_type), allocatable, intent(out) :: error
        type(PengRobinson) :: compound
        character(len=:), allocatable :: name
        real(wp) :: real_values(3), calc_values(3)

        name = "methane"

        compound = PR(name, ac=2.4885_wp, b=0.02664_wp, k=0.3923_wp)
        call compound%get_critical_constants()

        real_values = [191.1554, 46.4135, 0.0115]
        calc_values = [compound%tc, compound%pc, compound%w]

        call check(error, maxval(abs(real_values - calc_values)) > errmax)

        if (allocated(error)) return
    end subroutine test_get_critical

    subroutine test_attractive_parameter(error)
        type(error_type), allocatable, intent(out) :: error

        character(len=:), allocatable :: name
        type(PengRobinson) :: compounds(3)
        type(scalar_property) :: a_calc(3)
        real(wp) :: real_values(3)
        real(wp) :: temperature
        integer :: i

        real_values = [ &
            2.2159967041164599_wp, &
            6.6275760017729199_wp, &
            12.476490465688833_wp &
            ]
        a_calc = null_scalar_property(3)
        temperature = 250.0_wp

        compounds(1) = PR("methane", tc=191.15_wp, pc=46.41_wp, w=0.0115_wp)
        compounds(2) = PR("ethane", tc=305.3_wp, pc=49.0_wp, w=0.099_wp)
        compounds(3) = PR("propane", tc=369.9_wp, pc=42.5_wp, w=0.1521_wp)

        compounds%T = temperature
        a_calc = compounds%a_parameter()

        call check(error, maxval(abs(real_values - a_calc)) < ERRMAX)

        if (allocated(error)) return
    end subroutine test_attractive_parameter

    subroutine test_del1_definition(error)
        ! Use the correct delta 1 parameter
        type(error_type), allocatable, intent(out) :: error
        type(PengRobinson) :: compound

        compound = PR("name", 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0)
        call check(error, abs(1.0_wp + sqrt(2.0_wp) - compound%del1) < errmax)

        if (allocated(error)) return
    end subroutine test_del1_definition
end module
