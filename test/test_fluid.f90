module fixtures
    use feos
    implicit none

    private
    public PengRobinsonFluid

contains
    function PengRobinsonFluid() result(fluid)
        type(CubicFluid) :: fluid

        integer, parameter :: n = 3

        type(PengRobinson) :: compounds(n)
        type(ClassicVdW) :: mixing_rule
        type(binary_property) :: aij, bij

        real(wp) :: kij(n, n)
        real(wp) :: lij(n, n)

        compounds(1) = PR("methane", tc=191.15_wp, pc=46.41_wp, w=0.0115_wp)
        compounds(2) = PR("ethane",  tc=305.3_wp,  pc=49.0_wp,  w=0.099_wp)
        compounds(3) = PR("propane", tc=369.9_wp,  pc=42.5_wp,  w=0.1521_wp)

        kij = reshape(&
            [0.0, 0.4, 0.3, &
             0.4, 0.0, 0.1, &
             0.3, 0.1, 0.0],&
             [n, n] &
            )
        lij = 0*kij

        mixing_rule = ClassicVdW(kij=kij, lij=lij)

        compounds%moles = 1.0_wp/n

        fluid%components = compounds
        fluid%mixing_rule = mixing_rule
    end function
end module fixtures

module test_fluid
    use feos
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none

    private
    public :: collect_fluid

contains
    ! =============================================================================
    !  Tests that will be run
    ! -----------------------------------------------------------------------------
    subroutine collect_fluid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("FLUID: Helmholtz", test_ar), &
                    new_unittest("FLUID: vsolve", test_vsolve) &
                    ]
    end subroutine collect_fluid
    ! =============================================================================

    subroutine test_ar(error)
        use fixtures, only: PengRobinsonFluid
        type(error_type), allocatable, intent(out) :: error
        real(wp) :: v, t

        type(CubicFluid) :: fluid
        type(scalar_property) :: ar

        fluid = PengRobinsonFluid()
        
        v = 5_wp
        t = 150_wp
        ar = fluid%residual_helmholtz(v, t)   

        call check(error, abs(-1.2268112739280079 - ar%val) > 1.0d-8)
    end subroutine test_ar

    subroutine test_vsolve(error)
        use fixtures, only: PengRobinsonFluid
        type(error_type), allocatable, intent(out) :: error
        real(wp) :: v, t, p

        type(CubicFluid) :: fluid

        fluid = PengRobinsonFluid()

        p = 0.8
        v = fluid%vsolve(p, t)

        call check(error, abs(5.0348022921400115E-002 - v) > 1.0d-8)

    end subroutine test_vsolve
end module
