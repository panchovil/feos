module test_mixing
    use constants
    use properties
    use cubic_eos
    use mixing_rules
    use testdrive, only: new_unittest, unittest_type, error_type, check

    implicit none

contains

    subroutine test_classicvdw_def(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 2
        integer :: i, j
        type(ClassicVdW) :: vdw_mix

        character(len=:), allocatable :: name
        type(PengRobinson) :: compounds(n)
        type(scalar_property(n)) :: d, b, d1

        real(wp) :: kij(n, n)
        real(wp) :: lij(n, n)
        real(wp) :: T = 250.0_wp

        kij = reshape( &
              [0.0_wp, 0.4_wp, &
               0.4_wp, 0.0_wp], &
              [n, n] &
              )
        lij = reshape( &
              [0.0_wp, 0.6_wp, &
               0.6_wp, 0.0_wp], &
              [n, n] &
              )

        name = "methane"
        compounds(1) = PR("methane", ac=2.489_wp, b=0.0266_wp, k=0.3923_wp)
        compounds(1)%moles = 0.5_wp
        compounds(1)%T = T
        name = "ethane"
        compounds(2) = PR(name, ac=6.013_wp, b=0.0403_wp, k=0.5247_wp)
        compounds(2)%moles = 0.5_wp
        compounds(2)%T = T

        vdw_mix = ClassicVdW(kij=kij, lij=lij)

        do i = 1, 2
            print *, trim(compounds(i)%name), &
                compounds(i)%ac, compounds(i)%b, compounds(i)%k
        end do

        do i = 1, 2
            print *, vdw_mix%kij(:, i)
        end do
        call vdw_mix%mix(compounds, d, b, d1)
        print *, d%val, d%dn(:)

    end subroutine test_classicvdw_def
end module test_mixing

program main
    use test_mixing
    use testdrive, only: new_unittest, unittest_type, error_type, check
    type(error_type), allocatable :: error

    call test_classicvdw_def(error)
end program main
