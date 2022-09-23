module cubic_eos
   !! This module encompass the basic structure of equations of states.
   !! It considers the generic equation of 
   !! \[P = \frac{RT}{V-b} - \frac{a(T)}{()}\]

   use constants, only: wp, R
   use properties, only: scalar_property

   implicit none

   private

   ! Base objects
   public :: CubicEoS
   public :: PengRobinson
   public :: SoaveRedlichKwong
   public :: RKPR

   ! Objects factories
   public :: PR
   public :: SRK
   public :: RK_PR

   type :: CubicEoS
      !! Pure Compound parameters
      character(len=200) :: name    !! Compound name
      real(wp) :: T = 10000         !! Temperature
      real(wp) :: V = 20            !! Volume
      real(wp) :: ac                !! Critical atractive parameter
      real(wp) :: b                 !! Repulsive parameter
      real(wp) :: tc                !! Critical temperature
      real(wp) :: pc                !! Critical pressure
      real(wp) :: w                 !! Accentric factor
      real(wp) :: k                 !! Atractive parameter constant
      real(wp) :: del1              !! \[delta_1\] parameter
      real(wp) :: a = 0             !! Atractive parameter valuated at temperature
      real(wp) :: dadt = 0          !! Atractive parameter first derivative with tempetarue
      real(wp) :: da2dt2 = 0        !! Atractive parameter second derivative with tempetarue
   contains
      procedure :: a_parameter => a_classic
      procedure :: b_parameter => b_classic
   end type CubicEoS

   type, extends(CubicEoS) :: PengRobinson
       !! Peng-Robinson EoS
   contains
      procedure :: get_params => get_params_pr
      procedure :: get_critical_constants => get_critical_constants_pr
   end type

   type, extends(CubicEoS) :: SoaveRedlichKwong
       !! Soave-Redlich-Kwong EoS
   contains
      procedure :: get_params => get_params_srk
      procedure :: get_critical_constants => get_critical_constants_srk
   end type

   type, extends(CubicEoS) :: RKPR
       !! RKPR EoS
   contains
      procedure :: a_parameter => a_rkpr
      procedure :: get_params => get_params_rkpr
      procedure :: get_critical_constants => get_critical_constants_rkpr
   end type

contains
   ! ==========================================================================
   !  Peng Robinson EoS
   !  PR Equation of state methods
   ! --------------------------------------------------------------------------
   function PR(name, ac, b, tc, pc, w, k)
       !! PengRobinson factory
       type(PengRobinson) :: PR
       character(len=:), allocatable, intent(in) :: name
       real(wp),           optional, intent(in) :: ac
       real(wp),           optional, intent(in) :: b
       real(wp),           optional, intent(in) :: tc
       real(wp),           optional, intent(in) :: pc
       real(wp),           optional, intent(in) :: w
       real(wp),           optional, intent(in) :: k

       character(len=200) :: name_in
       real(wp) :: del1 = 1.0_wp + sqrt(2.0_wp)

       name_in = name

       if (present(ac) .and. present(b) .and. present(k)) then
           PR = PengRobinson(&
              name=name, ac=ac, b=b, tc=0, pc=0, w=0, k=k, del1=del1&
              )
           call PR%get_critical_constants()

       else if (present(tc) .and. present(pc) .and. present(w)) then
           PR = PengRobinson(&
              name=name, ac=0, b=0, tc=tc, pc=pc, w=w, k=0, del1=del1&
              )
           call PR%get_params()
       end if
   end function PR

   subroutine get_params_pr(self)
      class(PengRobinson) :: self
      real(wp) :: Zc, OMa, OMb, RT, Vceos
      RT = R*self%Tc
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%ac = OMa*RT**2/self%Pc
      self%b = OMb*RT/self%Pc
      Vceos = Zc*RT/self%Pc

      associate(k => self%k, w => self%w )
      ! k (or m) constant to calculate attractive parameter depending on temperature
      if (self%w <= 0.491) then
         ! m from PR
         k = 0.37464 + 1.54226*w - 0.26992*w**2
      end if
      if (w > 0.491) then
         ! PR78
         k = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
      end if
      end associate
   end subroutine get_params_pr
   
   subroutine get_critical_constants_pr(self)
      class(PengRobinson) :: self

      real(wp) :: OMa, OMb, Zc, al, be, ga, Vceos
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%Tc = OMb*self%ac/(OMa*R*self%b)
      self%Pc = OMb*R*self%Tc/self%b
      Vceos = Zc*R*self%Tc/self%Pc
         
      al = -0.26992
      be = 1.54226
      ga = 0.37464 - self%k
      self%w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
   end subroutine get_critical_constants_pr
   ! ==========================================================================

   ! ==========================================================================
   !  SoaveRedlichKwong EoS
   !  SoaveRedlichKwong methdos
   ! --------------------------------------------------------------------------
   function SRK(name, ac, b, tc, pc, w, k)
       !! SoaveRedlichKwong factory
       type(SoaveRedlichKwong) :: SRK
       character(len=200), optional, intent(in) :: name
       real(wp), optional, intent(in) :: ac
       real(wp), optional, intent(in) :: b
       real(wp), optional, intent(in) :: tc
       real(wp), optional, intent(in) :: pc
       real(wp), optional, intent(in) :: w
       real(wp), optional, intent(in) :: k

       real(wp) :: del1 = 1.0_wp

       if (present(ac) .and. present(b) .and. present(k)) then
           SRK = SoaveRedlichKwong(&
               name=name, ac=ac, b=b, tc=0, pc=0, w=0, k=k, del1=del1&
               )
           call SRK%get_critical_constants()
       else if (present(tc) .and. present(pc) .and. present(w)) then
           SRK = SoaveRedlichKwong(&
               name=name, ac=0, b=0, tc=tc, pc=pc, w=w, k=0, del1=del1&
               )
           call SRK%get_params()
       end if
   end function SRK

   subroutine get_params_srk(self)
      class(SoaveRedlichKwong) :: self
      
      real(wp) :: Zc, OMa, OMb, RT, Vceos
      rt = R*self%Tc
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%ac = OMa*rt**2/self%Pc
      self%b = OMb*rt/self%Pc
      Vceos = Zc*rt/self%Pc

      associate(k => self%k, w => self%w )
         ! k (or m) constant to calculate attractive parameter depending on temperature
         k = 0.48 + 1.574*w - 0.175*w**2
      end associate
   end subroutine get_params_srk

   subroutine get_critical_constants_srk(self)
      class(SoaveRedlichKwong) :: self

      real(wp) :: OMa, OMb, Zc, al, be, ga, Vceos
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%Tc = OMb*self%ac/(OMa*R*self%b)
      self%Pc = OMb*R*self%Tc/self%b
      Vceos = Zc*R*self%Tc/self%Pc
      al = -0.175
      be = 1.574
      ga = 0.48 - self%k
      self%w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
   end subroutine get_critical_constants_srk
   ! ==========================================================================
   
   ! ==========================================================================
   !  RKPR methods
   ! --------------------------------------------------------------------------
   function RK_PR(name, ac, b, tc, pc, w, k, del1)
       character(len=200), optional, intent(in) :: name
       real(wp), optional, intent(in) :: ac
       real(wp), optional, intent(in) :: b
       real(wp), optional, intent(in) :: tc
       real(wp), optional, intent(in) :: pc
       real(wp), optional, intent(in) :: w
       real(wp), optional, intent(in) :: k
       real(wp), optional, intent(in) :: del1

       type(RKPR) :: RK_PR

       if (present(ac) .and. present(b) .and. present(k)) then
           RK_PR = RKPR(name=name, ac=ac, b=b, tc=0, pc=0, w=0, k=k, del1=del1)
           call RK_PR%get_critical_constants()
       else if (present(tc) .and. present(pc) .and. present(w)) then
           RK_PR = RKPR(name=name, ac=0, b=0, tc=tc, pc=pc, w=w, k=0, del1=del1)
           call RK_PR%get_params()
       end if
   end function RK_PR

   subroutine get_params_rkpr(self)
      class(RKPR) :: self
      real(wp) :: Zc, OMa, OMb, RT, Vceos
      RT = R*self%Tc
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%ac = OMa*RT**2/self%Pc
      self%b = OMb*RT/self%Pc
      Vceos = Zc*RT/self%Pc
      print *, "EXCEPTION: Not implemented error"
   end subroutine get_params_rkpr
   
   subroutine get_critical_constants_rkpr(self)
      class(RKPR) :: self
      real(wp) :: Zc, OMa, OMb, RT, Vceos
      RT = R*self%Tc
      call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

      self%ac = OMa*RT**2/self%Pc
      self%b = OMb*RT/self%Pc
      Vceos = Zc*RT/self%Pc
      print *, "EXCEPTION: Not implemented error"
   end subroutine get_critical_constants_rkpr
   ! ==========================================================================

   ! ==========================================================================
   !  Attractive terms subroutines
   ! --------------------------------------------------------------------------
   subroutine a_classic(self)
      !! Calculate the atractive parameter at T temperature.
      class(CubicEoS) :: self

      real(8) :: Tr
      real(8) :: ac

     associate(T => self%T, Tc => self%Tc, k  => self%k)
         Tr = T/Tc
         self%a = ac*(1 + k*(1 - sqrt(Tr)))**2
         self%dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
         self%da2dT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5)
      end associate
   end subroutine a_classic
   
   subroutine a_rkpr(self)
      !! Calculate the atractive parameter at T temperature.
      class(RKPR) :: self

      real(8) :: Tr
      real(8) :: ac
      
      associate(T => self%T, Tc => self%Tc, k  => self%k)
         Tr = T/Tc
         self%a = ac*(3/(2 + Tr))**k
         self%dadT = -k*self%a/Tc/(2 + Tr)
         self%da2dT2 = -(k + 1)*self%dadT/Tc/(2 + Tr)
      end associate
   end subroutine a_rkpr
   ! ==========================================================================
   
   ! ==========================================================================
   !  Repulsive term routines
   ! --------------------------------------------------------------------------
   subroutine b_classic(self)
      !! Classic CubicEoS has a constant repulsive parameter
      class(CubicEoS) :: self
      self%b = self%b
   end subroutine b_classic
   ! ==========================================================================

   ! ==========================================================================
   !  Auxiliar routines
   ! --------------------------------------------------------------------------
   subroutine get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
      !! Calculate Zc, OMa and OMb from the delta_1 parameter.
      real(wp), intent(in)  :: del1 !! delta_1 parameter
      real(wp), intent(out) :: Zc   !! Critical compressibility factor
      real(wp), intent(out) :: OMa  !! OMa
      real(wp), intent(out) :: OMb  !! OMb

      real(wp) :: d1, y

      d1 = (1.d0 + del1**2.d0)/(1.d0 + del1)
      y = 1.d0 + (2.d0*(1.d0 + del1))**(1.0d0/3.d0) + (4.d0/(1.d0 + del1))**(1.0d0/3)
      OMa = (3.d0*y*y + 3.d0*y*d1 + d1**2.d0 + d1 - 1.0d0)/(3.d0*y + d1 - 1.0d0)**2.d0
      OMb = 1.d0/(3.d0*y + d1 - 1.0d0)
      Zc = y/(3.d0*y + d1 - 1.0d0)
   end subroutine get_Zc_OMa_OMb
   ! ==========================================================================
   
   subroutine residual_helmholtz(nc, V, T, D, D1, Bmix, Ar)
      !! Mixture of fluids helmholtz energy
      integer, intent(in) :: nc !! Number of components
      real(wp), intent(in) :: V !! Volume
      real(wp), intent(in) :: T !! Temperature

      real(wp), intent(in) :: D !! Atractive parameter times moles
      real(wp), intent(in) :: D1 !! Delta_1 parameter
      real(wp), intent(in) :: Bmix !! Repulsive parameter
      type(scalar_property(nc)), intent(out) :: Ar !! Residual Helmholtz energy object

      integer :: i, j

      real(wp) :: totn, D2
      real(wp) :: f, g, fv, gv, fB, fv2, gv2

      real(wp) :: auxD2, fD1, fBD1, fVD1, fD1D1
      real(wp) :: AUX, FFB, FFBV, FFBB

      real(wp) :: dD1i(nc), dD1ij(nc, nc)

      real(wp) :: dDdT, dDdT2, dDi(nc), dDiT(nc), dDij(nc, nc)

      real(wp) :: dBi(nc), dBij(nc, nc)

      D2 = (1._wp - D1)/(1._wp + D1)

      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = R*log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = R*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

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
      Ar%val = -TOTN*g*T - D*f
      Ar%dv = -TOTN*gv*T - D*fv
      Ar%dv2 = -TOTN*gv2*T - D*fv2

      AUX = R*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      do i = 1, nc
         Ar%dn(i) = -g*T + FFB*dBi(i) - f*dDi(i) - D*fD1*dD1i(i)
         Ar%dvn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i) - D*fVD1*dD1i(i)
         do j = 1, i
            Ar%dn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                         + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
            Ar%dn2(i, j) = Ar%dn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                         - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                         - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
            Ar%dn2(j, i) = Ar%dn2(i, j)
         end do
      end do

      ! TEMPERATURE DERIVATIVES
      Ar%dt = -TOTN*g - dDdT*f
      Ar%dtv = -TOTN*gv - dDdT*fV
      Ar%dt2 = -dDdT2*f
      do i = 1, nc
         Ar%dtn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i) - dDdT*fD1*dD1i(i)
      end do
   end subroutine residual_helmholtz
end module cubic_eos
