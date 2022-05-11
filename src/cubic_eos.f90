module cubic_eos
    !! Module that encompass the calculations of the residual Helmholtz energy
    !! and related properties like fugacity coefficents.
    use constants
    implicit none

    contains

    subroutine kij_tdep(T, kij, dkijdt, dkij2dt2)
        !! Kij with temperature dependance according to the equation:
        !! \[ K_{ij}(T) = K_{ij\infty} + K_{ij0} e^{T/T^*} \]
        !! The parameters of the equation are obtained from the mixture module
        use mixture, only: nc, kij_0, kij_inf, T_star
        implicit none
        real(wp), intent(in) :: T !! Temperature
        real(wp), intent(out) :: kij(nc, nc) !! Binary interaction parameter matrix
        real(wp), intent(out) :: dkijdt(nc, nc) !! Binary interaction parameter first derivative with T matrix
        real(wp), intent(out) :: dkij2dt2(nc, nc) !! Binary interaction parameter second derivative with T matrix

        integer :: i

        do i = 1, nc
            kij(:i - 1, i) = kij_inf(:i - 1, i) + kij_0(:i - 1, i)*exp(-T/T_star(:i - 1, i))
            dkijdt(:i - 1, i) = -kij_0(:i - 1, i)/T_star(:i-1,i)*exp(-T/T_star(:i-1,i))
            dkij2dt2(:i - 1, i) = kij_0(:i - 1, i)/T_star(:i-1,i)**2*exp(-T/T_star(:i-1,i))
        end do
    end subroutine kij_tdep

    subroutine aTder(T, a, dadT, dadT2)
        !! Calculate the atractive parameter at T temperature.
        !! the subroutine will read the mixture's model and based on that
        !! will use the corresponding rule.
        use mixture, only: nmodel, ac, k, tc, nc
        real(wp), intent(in) :: T  !! Temperature

        real(wp), intent(out) :: a(nc)  !! Atractive parameter at T
        real(wp), intent(out) :: dadT(nc) !! First deritvative with T
        real(wp), intent(out) :: dadT2(nc)  !! Second derivative with T

        real(wp) :: Tr(nc) ! Reduced temperature

        Tr = T/Tc
        if (nmodel <= 2) then
            ! SRK and PR
            a = ac*(1 + k*(1 - sqrt(Tr)))**2
            dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
            dadT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5)
        else if (nmodel == 3) then
            ! RKPR EOS
            a = ac*(3/(2 + Tr))**k
            dadT = -k*a/Tc/(2 + Tr)
            dadT2 = -(k + 1)*dadT/Tc/(2 + Tr)
        end if
    end subroutine aTder

    subroutine aijTder(T, aij, daijdT, daij2dT2)
        !! Calculate the binary atractive term matrix
        use mixture, only: nc, ntdep, kij, kij_0, lij, ncomb
        use mixing_rules, only: quadratic
        implicit none

        real(wp), intent(in) :: T !! Temperature

        real(wp), intent(out) :: aij(nc, nc) !! Atractive binary terms matrix
        real(wp), intent(out) :: daijdT(nc, nc) !! Atractive binary terms matrix first derivative with temperature
        real(wp), intent(out) :: daij2dT2(nc, nc) !! Atractive binary terms matrix second derivative with temperature

        real(wp) :: dkijdt(nc, nc), dkij2dt2(nc, nc) ! kij T derivatives
        real(wp) :: a(nc), dadT(nc), da2dT2(nc) ! Atractive parameter and T derivatives
        real(wp) :: b(nc), bij(nc, nc) ! Repulsive parameter (just to use as input in subroutine)

        b = 0  ! Here only the aij for the mixture, so there is no need to use the real b

        select case (ntdep)
            case (1)
                ! Kij exponential temperature dependance
                call kij_tdep(T, kij, dkijdt, dkij2dt2)
            case default
                kij = kij_0
                dkijdt = 0
                dkij2dt2 = 0
        end select 

        ! Calculate pure compound atractive parameters at T
        call aTder(T, a, dadT, da2dT2)

        ! Apply combining rule
        select case (ncomb)
            case default
                call quadtratic(a, b, kij, dadt, da2dt2,  dkijdt, dkij2dt2, &
                    lij, aij, daijdt, daij2dt2, bij, nc)
        end select
    end subroutine aijTder

    subroutine aTnder(T, n, a, dadni, da2dniT2, da2dnij2, dadT, da2dT2)
       use mixture, only: nc
       implicit none

       real(wp), intent(in) :: T !! Temperature
       real(wp), intent(in) :: n(nc) !! Matrix with number of moles

       real(wp), intent(out) :: a !! Mixture atractive parameter
       real(wp), intent(out) :: dadni(nc) !! Atractive parameter first derivative with number of moles
       real(wp), intent(out) :: da2dniT2(nc) !! Atractive parameter second derivative with moles and Tempeature
       real(wp), intent(out) :: da2dnij2(nc, nc) !! Atractive parameter second derivative with number of moles

       real(wp), intent(out) :: dadT(nc) !! Atractive parameter first derivative with Tempeature
       real(wp), intent(out) :: da2dT2(nc) !! Atractive parameter first derivative with Tempeature

       real(wp) :: aij(nc, nc), daijdT(nc, nc), daijdT2(nc, nc)
       real(wp) :: aux, aux2
       integer :: i, j

       call aijTder(T, aij, daijdT, daijdT2)
       a = 0.0_wp
       dadT = 0.0_wp
       da2dT2 = 0.0_wp
       do i = 1, nc
          aux = 0.0_wp
          aux2 = 0.0_wp
          dadni(i) = 0.0_wp
          da2dniT2(i) = 0.0_wp
          do j = 1, nc
             dadni(i) = dadni(i) + 2*n(j)*aij(i, j)
             da2dniT2(i) = da2dniT2(i) + 2*n(j)*daijdT(i, j)
             da2dnij2(i, j) = 2*aij(i, j)

             aux = aux + n(j)*aij(i, j)
             aux2 = aux2 + n(j)*daijdT2(i, j)
          end do
          a = a + n(i)*aux
          dadT = dadT + n(i)*da2dniT2(i)/2
          da2dT2 = da2dT2 + n(i)*aux2
       end do
    end subroutine aTnder
end module cubic_eos
