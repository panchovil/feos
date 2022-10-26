module cubic_eos_core
    !! Abstract interfaces for Cubic Equations of State.
   implicit none

   private
   public :: BaseCEOS

   type, abstract :: BaseCEOS
       !! Abstract class that represents a basic Cubic Equation Of State.
       !!
       !! This derived type assumes a Cubic Equation of State uses 
       !! two parameters that can be calculated with the attributes 
       !! contained inside the derived type. It only enforces that there *must*
       !! be a method that calculates the attractive and repulsive parameters.
       !! Adding extra methods for the calculation of extra parameters won't
       !! be a problem, but it must be taken care that the mixing rule and
       !! helmholtz calculation functions take that into account.
   contains
      !> Attractive parameter calculator
      procedure(a_parameter), deferred :: a_parameter 
      !> Repulsive parameter calculator
      procedure(b_parameter), deferred :: b_parameter
   end type BaseCEOS

   abstract interface
      elemental function a_parameter(self) result(a)
         use properties, only: scalar_property
         import BaseCEOS
         class(BaseCEOS), intent(in) :: self !! Cubic EoS derived type.
         type(scalar_property) :: a !! Calculated a parameter.
      end function
      elemental function b_parameter(self) result(b)
         use properties, only: scalar_property
         import BaseCEOS
         class(BaseCEOS), intent(in) :: self !! Cubic EoS derived type.
         type(scalar_property) :: b !! Calculated output parameter.
      end function
   end interface
end module cubic_eos_core
