module fluid
    !! Fluid module
    !! 
    !! This module contains the derived type that represents a fluid. It's
    !! considered as a mixture of some kind of equation of state 
    !! (pure components), some mixing rule that it uses, and a method
    !! to calculate the residual helmholtz energy.
    use constants
    use properties
    use cubic_eos
    use mixing_rules
    use fluid_core

    implicit none
    private

    public :: CubicFluid

    type, extends(BaseFluid) :: CubicFluid
        !! Cubic fluid.
        !!
        !! Derived type that represents a fluid that uses a Cubic EOS.
        !! It contains an array of pure components (that are represented by
        !! a Cubic EOS) and the corresponding mixing rule that's used.
        real(wp) :: p
        real(wp) :: v
        real(wp) :: t
        class(CubicEoS), allocatable :: components(:) !! Components
        class(ClassicVdW), allocatable :: mixing_rule !! Mixing rule
    contains
        !> Residual Helmholtz energy
        procedure :: residual_helmholtz => residual_helm
        !> Change fluid's PVT state
        procedure :: set_to => set_to 
        !> Copy the fluid
        procedure :: copy => copy 
    end type CubicFluid

contains

    pure subroutine set_to(self, pressure, volume,  temperature)
        !! Change a fluid PVT state.
        !!
        !! Works based on specification. Iterative procedures will be
        !! called in case of a PV or PT specification is made.

        class(CubicFluid), intent(in out) :: self
        real(wp), optional, intent(in) :: pressure !! Pressure [bar]
        real(wp), optional, intent(in) :: volume !! Volume [dm^2]
        real(wp), optional, intent(in) :: temperature !! Temperature [K]

        if (present(volume)) then
            self%components%v = volume
        end if

        if (present(temperature)) then
            self%components%t = temperature
        end if
    end subroutine

    pure function copy(self) result(new)
        !! Copy a fluid.
        !!
        !! Method that returns a copy of a fluid.
        class(CubicFluid), intent(in) :: self !! Fluid to copy
        type(CubicFluid), allocatable :: new !! Fluid copy

        allocate(CubicFluid :: new)
        new%components = self%components
        new%mixing_rule = self%mixing_rule
    end function
    
    elemental function residual_helm(self, v, t) result(ar)
        !! Residual Helmholtz Energy.
        !!
        !! This method calculates the residual Helholtz energy and it's relevant
        !! derivatives.
        !!
        !! Internally defines a temporary fluid at the desired temperature
        !! and volume values, call it's defined mixing rule and calculate
        !! the residual Helmholtz energy with it's derivatives with the
        !! generic cubic residual Helmoltz function.
        !!
        !! The residual Helmholtz energy is calculated with the RKPR EOS [1]
        !! approach, which uses a a third parameter (delta 1)
        !!\[ \alpha^r(T, V, n) = -n \ln(1 - B/V) - \frac{D(T)}{RTB(\delta_1 - \delta_2)}\ln\left(\frac{1+\delta_1 B/V}{1+\delta_2 B/V}\right) \]
        !!
        !! with \[\delta_2 = \frac{1 - \delta_1}{1 + \delta_2}\]

        class(CubicFluid), intent(in) :: self
        !> Volume [dm^3]
        real(wp), intent(in) :: v 
        !> Temperature [K]
        real(wp), intent(in) :: t 
        !> Residual helmholtz energy and it's derivatives
        type(scalar_property) :: ar 

        ! Internal variables
        ! concentrations vector
        real(wp) :: moles(size(self%components))
        ! Set of compounds
        class(CubicEoS), allocatable :: compounds(size(self%components))
        ! Fluid copy, this is used to avoid modifications of the original
        ! fluid.
        class(CubicFluid), allocatable :: fluid_tmp
        ! Number of components
        integer :: n

        ! Mixture "as a fluid" parameters.
        type(scalar_property) :: d, b, d1

        n = size(self%components)
        ar = null_scalar_property(n)
        fluid_tmp = self%copy()

        ! Go to desired volume and temperature
        call fluid_tmp%set_to(volume=v, temperature=t)
    
        ! Copy compounds into internal variables 
        compounds = fluid_tmp%components
        moles = fluid_tmp%components%moles

        ! Call fluid mixing rule
        call fluid_tmp%mixing_rule%mix(compounds, d, b, d1)
        
        ar = rkpr_residual_helmholtz(n, moles, v, t, d, d1, b)
    end function residual_helm

end module fluid
