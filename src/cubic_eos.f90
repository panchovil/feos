module cubic_eos
   !! Module that encompass the calculations of the residual Helmholtz energy
   !! and related properties like fugacity coefficents.

   use constants

   implicit none

   private
   public :: pure_compound
   public :: pr
   public :: srk
   public :: rkpr

   type, abstract :: pure_compound
      !! Generic EoS
      character(len=:), allocatable :: name !! Compound name
      real(wp) :: ac !! Critical atractive parameter
      real(wp) :: b !! Repulsive parameter
      real(wp) :: tc !! Critical temperature
      real(wp) :: pc ! Critical pressure
      real(wp) :: w !! Accentric factor
      real(wp) :: k !! Atractive parameter constant
      real(wp) :: a = 0 !! Atractive parameter valuated at temperature
      real(wp) :: dadt = 0 !! Atractive parameter first derivative with tempetarue
      real(wp) :: da2dt2 = 0 !! Atractive parameter second derivative with tempetarue
   contains
      procedure :: a_t => a_parameter
      procedure :: get_params => get_params
      procedure :: get_critical_constants => get_critical_constants
   end type pure_compound

   type, extends(pure_compound) :: pr
      !! Peng-Robinson EoS
      real(wp) :: del1 = 1.0_wp + sqrt(2.0_wp) !! \[\delta_1\] parameter.
   end type pr

   type, extends(pure_compound) :: srk
      !! Soave-Redlich-Kwong EoS
      real(wp) :: del1 = 1 !! \[\delta_1\] parameter.
   end type srk

   type, extends(pure_compound) :: rkpr
      real(wp) :: del1
   end type rkpr

contains

   subroutine get_params(self)
      class(pure_compound) :: self

      real(wp) :: Zc, OMa, OMb, RT, Vceos

      RT = R*self%Tc

      associate( k => self%k, w => self%w )
      select type (self)
      class is (pr)
         call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)
         self%ac = OMa*RT**2/self%Pc
         self%b = OMb*RT/self%Pc
         Vceos = Zc*RT/self%Pc
         ! m constant to calculate a depending on temperature
         if (w <= 0.491) then
            ! m from PR
            k = 0.37464 + 1.54226*w - 0.26992*w**2
         end if
         if (w > 0.491) then
            ! PR78
            k = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
         end if

      class is (srk)
         call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)
         self%ac = OMa*RT**2/self%Pc
         self%b = OMb*RT/self%Pc
         Vceos = Zc*RT/self%Pc
         k = 0.48 + 1.574*w - 0.175*w**2

      class is (rkpr)
         call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)
         self%ac = OMa*RT**2/self%Pc
         self%b = OMb*RT/self%Pc
         Vceos = Zc*RT/self%Pc
         print *, "Not implemented yet"
      end select
      end associate
   end subroutine get_params


   subroutine get_critical_constants(self)
      class(pure_compound) :: self
      
      real(wp) :: OMa, OMb, Zc, al, be, ga, Vceos

      select type (self)
      class is (pr)
         call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

         self%Tc = OMb*self%ac/(OMa*R*self%b)
         self%Pc = OMb*R*self%Tc/self%b
         Vceos = Zc*R*self%Tc/self%Pc

         al = -0.26992
         be = 1.54226
         ga = 0.37464 - self%k

         self%w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
      class is (srk)
         call get_Zc_OMa_OMb(self%del1, Zc, OMa, OMb)

         self%Tc = OMb*self%ac/(OMa*R*self%b)
         self%Pc = OMb*R*self%Tc/self%b
         Vceos = Zc*R*self%Tc/self%Pc

         al = -0.175
         be = 1.574
         ga = 0.48 - self%k
         self%w = 0.5*(-be + sqrt(be**2 - 4*al*ga))/al
      class is (rkpr)
         print *, "Not implemented yet"
      end select
   end subroutine get_critical_constants

   
   subroutine get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
      real(wp), intent(in) :: del1 !! RKPR delta_1 parameter
      real(wp), intent(out) :: Zc !! Critical compressibility factor
      real(wp), intent(out) :: OMa !! 
      real(wp), intent(out) :: OMb !!

      real(wp) :: d1, y

      d1 = (1.d0 + del1**2.d0)/(1.d0 + del1)
      y = 1.d0 + (2.d0*(1.d0 + del1))**(1.0d0/3.d0) + (4.d0/(1.d0 + del1))**(1.0d0/3)
      OMa = (3.d0*y*y + 3.d0*y*d1 + d1**2.d0 + d1 - 1.0d0)/(3.d0*y + d1 - 1.0d0)**2.d0
      OMb = 1.d0/(3.d0*y + d1 - 1.0d0)
      Zc = y/(3.d0*y + d1 - 1.0d0)
   end subroutine get_Zc_OMa_OMb


   subroutine a_parameter(self, T)
      !! Calculate the atractive parameter at T temperature.
      !! the subroutine will read the mixture's model and based on that
      !! will use the corresponding rule.
      class(pure_compound) :: self
      real(8), intent(in) :: T !! Temperature where to calculate

      real(8) :: Tr
      real(8) :: ac, k, Tc

      Tc = self%Tc
      ac = self%ac
      k = self%k
      Tr = T/Tc

      select type (self)
      class is (rkpr)
         self%a = ac*(3/(2 + Tr))**k
         self%dadT = -k*self%a/Tc/(2 + Tr)
         self%da2dT2 = -(k + 1)*self%dadT/Tc/(2 + Tr)

      class default
         self%a = ac*(1 + k*(1 - sqrt(Tr)))**2
         self%dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
         self%da2dT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5)
      end select
   end subroutine a_parameter

end module cubic_eos
