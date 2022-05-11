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
         dkijdt(:i - 1, i) = -kij_0(:i - 1, i)/T_star(:i - 1, i)*exp(-T/T_star(:i - 1, i))
         dkij2dt2(:i - 1, i) = kij_0(:i - 1, i)/T_star(:i - 1, i)**2*exp(-T/T_star(:i - 1, i))
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
         call quadtratic(nc, a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, &
                         lij, aij, daijdt, daij2dt2, bij)
      end select
   end subroutine aijTder

   subroutine amixTnder(nc, T, n, a, dadni, da2dniT2, da2dnij2, dadT, da2dT2)
      !! Atractive parameter of a mixture and it's derivatives
      integer, intent(in) :: nc !! Number of components
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

      ! TODO: An already calculated aij matrix could be the input
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
   end subroutine amixTnder

   subroutine d1nder(nc, n, d1, d1_mix, dD1dni, dD12dnij2)
      !! Delta_1 of the mixture and compositional derivatives.
      integer, intent(in) :: nc !! Number of components.
      real(wp), intent(in) :: n(nc) !! Array of mole numbers for each component.
      real(wp), intent(in) :: d1(nc) !! Array of delta 1 parameters for each component.
      real(wp), intent(out) :: d1_mix !! ??
      real(wp), intent(out) :: dD1dni(nc) !! delta1 parameter first derivative with composition.
      real(wp), intent(out) :: dD12dnij2(nc, nc) !! delta1 parameter second derivative with composition.

      integer :: i, j
      real(wp) :: totn ! Total number of moles

      d1_mix = 0.0_wp

      do i = 1, nc
         d1_mix = d1_mix + n(i)*d1(i)
      end do
      totn = sum(n)
      d1_mix = d1_mix/totn
      do i = 1, nc
         dD1dni(i) = (d1(i) - d1_mix)/totn
         do j = 1, nc
            dD12dnij2(i, j) = (2.0_wp*d1_mix - d1(i) - d1(j))/totn**2
         end do
      end do
   end subroutine d1nder

   subroutine Bmixnder(nc, n, bij, B_mix, dBdni, dB2dnij2)
      !! Repulsive parameter of the mixture and it's compositional derivatives.
      integer, intent(in) :: nc
      real(wp), intent(in) :: n(nc)
      real(wp), intent(in) :: bij(nc, nc)

      real(wp), intent(out) :: B_mix
      real(wp), intent(out) :: dBdni(nc)
      real(wp), intent(out) :: dB2dnij2(nc, nc)

      real(wp) :: totn, aux(nc)
      integer :: i, j

      totn = sum(n)
      B_mix = 0.0_wp
      aux = 0.0_wp

      do i = 1, nc
         do j = 1, nc
            aux(i) = aux(i) + n(j)*bij(i, j)
         end do
         B_mix = B_mix + n(i)*aux(i)
      end do

      B_mix = B_mix/totn

      do i = 1, nc
         dBdni(i) = (2*aux(i) - B_mix)/totn
         do j = 1, i
            dB2dnij2(i, j) = (2*bij(i, j) - dBdni(i) - dBdni(j))/totn
            dB2dnij2(j, i) = dB2dnij2(i, j)
         end do
      end do
   end subroutine Bnder

   subroutine helmholtz_energy(nc, ND, NT, n, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: ND, NT
      real(wp), intent(in) :: n(nc) !! Number of moles
      real(wp), intent(in) :: V !! 
      real(wp), intent(in) :: T

      real(wp), intent(out) :: Ar
      real(wp), intent(out) :: ArV
      real(wp), intent(out) :: ArTV
      real(wp), intent(out) :: ArV2
      real(wp), intent(out) :: Arn(nc)
      real(wp), intent(out) :: ArVn(nc)
      real(wp), intent(out) :: ArTn(nc)
      real(wp), intent(out) :: Arn2(nc, nc)

      real(wp) :: dBi(nc), dBij(nc, nc)
      real(wp) :: dDi(nc), dDij(nc, nc), dDiT(nc)
      real(wp) :: aij(nc, nc), daijdT(nc, nc), daijdT2(nc, nc)
      real(wp) :: Kij(nco, nco)
      real(wp) :: ac(nco), b(nco), del1(nco), rm(nco)

      TOTN = sum(rn)
      D1 = del1(1)
      D2 = (1 - D1)/(1 + D1)

      if (ncomb .lt. 2) then
         call amixTnder(nc, T, n, a, dadni, da2dniT2, da2dnij2, dadT, da2dT2)
         call Bmixnder(nc, rn, Bmix, dBi, dBij)
      else
         ! call Bcubicnder(nc,rn,Bmix,dBi,dBij)
         ! call DCubicandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      end if

      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = RGAS*log(1 - Bmix/V)

      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      gv = RGAS*Bmix/(V*(V - Bmix))

      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = RGAS*(1/V**2 - 1/(V - Bmix)**2)

      fB = -(f + V*fv)/Bmix

      ! Reduced Helmholtz Energy and derivatives
      Ar = -TOTN*g*T - D*f
      ArV = -TOTN*gv*T - D*fv
      ArV2 = -TOTN*gv2*T - D*fv2

      AUX = RGAS*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      do i = 1, nc
         Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i)
         ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i)
         if (ND .EQ. 2) then
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                            + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end if
      end do

      ! TEMPERATURE DERIVATIVES
      IF (NT .EQ. 1) THEN
         ArT = -TOTN*g - dDdT*f
         ArTV = -TOTN*gv - dDdT*fV
         ArTT = -dDdT2*f
         do i = 1, nc
            ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i)
         end do
      END IF
   end subroutine helmholtz_energy

   subroutine helmholtz_energy_RKPR(nco, NDE, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! Calculate the reduced residual Helmholtz Energy and it's derivatives with the RKPR EOS 
      IMPLICIT DOUBLE PRECISION(A - H, O - Z)
      PARAMETER(RGAS=0.08314472d0)
      dimension :: rn(nco), Arn(nco), ArVn(nco), ArTn(nco), Arn2(nco, nco)
      dimension dBi(nco), dBij(nco, nco), dD1i(nco), dD1ij(nco, nco)
      dimension dDi(nco), dDij(nco, nco), dDiT(nco)
      dimension aij(nco, nco), daijdT(nco, nco), daijdT2(nco, nco)
      COMMON/rule/ncomb

      nc = nco
      TOTN = sum(rn)
      call DELTAnder(nc, rn, D1, dD1i, dD1ij)
      D2 = (1 - D1)/(1 + D1)

      if (ncomb .lt. 2) then
         call Bnder(nc, rn, Bmix, dBi, dBij)
         call DandTnder(NTD, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
      else
         ! call Bcubicnder(nc,rn,Bmix,dBi,dBij)
         ! call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      end if

      !  The f's and g's used here are for Ar, not F (reduced Ar)
      !  This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = RGAS*log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = RGAS*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = RGAS*(1/V**2 - 1/(V - Bmix)**2)

      ! DERIVATIVES OF f WITH RESPECT TO DELTA1
      auxD2 = (1 + 2/(1 + D1)**2)
      fD1 = (1/(V + D1*Bmix) + 2/(V + D2*Bmix)/(1 + D1)**2) - f*auxD2
      fD1 = fD1/(D1 - D2)
      fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/(1 + D1)**2)
      fBD1 = fBD1/(D1 - D2)
      fVD1 = -(fV*auxD2 + 1/(V + D1*Bmix)**2 + 2/(V + D2*Bmix)**2/(1 + D1)**2)/(D1 - D2)
      fD1D1 = 4*(f - 1/(V + D2*Bmix))/(1 + D1)**3 + Bmix*(-1/(V + D1*Bmix)**2  & 
              + 4/(V + D2*Bmix)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
      fD1D1 = fD1D1/(D1 - D2)

      ! Reduced Helmholtz Energy and derivatives
      Ar = -TOTN*g*T - D*f
      ArV = -TOTN*gv*T - D*fv
      ArV2 = -TOTN*gv2*T - D*fv2

      AUX = RGAS*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      do i = 1, nc
         Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i) - D*fD1*dD1i(i)
         ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i) - D*fVD1*dD1i(i)
         if (NDE .EQ. 2) then
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                            + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                            - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                            - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end if
      end do

      ! TEMPERATURE DERIVATIVES
      IF (NTD .EQ. 1) THEN
         ArT = -TOTN*g - dDdT*f
         ArTV = -TOTN*gv - dDdT*fV
         ArTT = -dDdT2*f
         do i = 1, nc
            ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i) - dDdT*fD1*dD1i(i)
         end do
      END IF
   end subroutine helmholtz_energy_RKPR
end module cubic_eos
