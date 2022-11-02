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
        !> Pressure calculation
        procedure :: pressure => pressure
        !> lnphi calculation
        procedure :: lnphi => fugacity
        !> Volume solver
        procedure :: vsolve => solve_volume
    end type CubicFluid

contains

    pure subroutine set_to(self, pressure, volume, temperature)
        !! Change a fluid PVT state.
        !!
        !! Works based on specification. Iterative procedures will be
        !! called in case of a PV or PT specification is made.

        class(CubicFluid), intent(in out) :: self
        real(wp), optional, intent(in) :: pressure !! Pressure [bar]
        real(wp), optional, intent(in) :: volume !! Volume [dm^2]
        real(wp), optional, intent(in) :: temperature !! Temperature [K]

        if (present(volume)) then
            self%v = volume
            self%components%v = volume
        end if

        if (present(temperature)) then
            self%t = temperature
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
        new%p = self%p
        new%v = self%v
        new%t = self%t
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
        real(wp), allocatable :: moles(:)
        ! Set of compounds
        class(CubicEoS), allocatable :: compounds(:)
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

    elemental function pressure(self, v, t) result(p)
        !! Pressure calculation
        class(CubicFluid), intent(in) :: self
        !> Volume [dm^3]
        real(wp), intent(in) :: v 
        !> Temperature [K]
        real(wp), intent(in) :: t 
        !> Pressure and it's derivatives
        type(scalar_property) :: p 
        
        type(scalar_property) :: ar 

        integer :: n

        real(wp) :: moles
        
        n = size(self%components)
        moles = sum(self%components%moles)

        p = null_scalar_property(n)

        ar = self%residual_helmholtz(v, t)

        p%val = -ar%dv + moles * R * t / v
        p%dt = -ar%dtv + moles * R / v
        p%dv = -ar%dv2 - moles * R * t / v**2
        p%dn = -ar%dvn + R * t / v
    end function pressure
    
    pure function fugacity(self, v, t) result(lnphi)
        !! Pressure calculation
        class(CubicFluid), intent(in) :: self
        !> Volume [dm^3]
        real(wp), intent(in) :: v 
        !> Temperature [K]
        real(wp), intent(in) :: t 
        !> Pressure and it's derivatives
        type(scalar_property), allocatable :: lnphi(:)
        
        type(scalar_property) :: ar, p

        integer :: n

        real(wp) :: moles, RT

        n = size(self%components)
        moles = sum(self%components%moles)

        RT = R * t
        ar = null_scalar_property(n)
        p = null_scalar_property(n)
        allocate(lnphi(n))
        lnphi = null_scalar_property(n)

        ar = self%residual_helmholtz(v, t)
        p = self%pressure(v, t)

        lnphi%val = ar%dn/RT - log((p*v)/(RT))
        lnphi%dt = ar%dtn/RT
    end function fugacity

    recursive pure function solve_volume(self, p_obj, t, root, max_it) result(v)
        !! Solve volume root for a specified pressure at the fluid's temperature

        !> Fluid
        class(CubicFluid), intent(in) :: self
        !> Objective pressure
        real(wp), intent(in) :: p_obj
        !> Temperature
        real(wp), intent(in) :: t
        !> Desired root, could be either "vapor" or "liquid"
        character(len=*), optional, intent(in) :: root 
        !> Max number of iterations, defaults to 50
        integer, optional, intent(in) :: max_it
        !> Output volume
        real(wp) :: v
        !> fluid at new volume to use inside iterations
        class(CubicFluid), allocatable :: fluid

        ! Mixture parameters and calculated pressure
        type(scalar_property) :: d, b, d1, p
        
        ! Inside variables
        real(wp) :: error, v0, z_min, z_max, z
        integer :: it, maxit

        if (.not. present(max_it)) then
             maxit = 50
        else
            maxit = max_it
        end if
        
        fluid = self%copy()

        call fluid%mixing_rule%mix(fluid%components, d, b, d1)

        z_min = 0
        z_max = 1 - 0.01 * t / (10000*b)

        if (present(root)) then
            select case (root)
            case ("vapor")
                z = min(0.5_wp, b * p_obj / (r * t))
                v0 = b/z
            case ("liquid")
                z = 0.5_wp
                v0 = b/z
            end select
            p = fluid%pressure(v0, t)
            error = abs(p - p_obj)

            it = 0
            do while(error >= errmax .and. it < maxit)
                it = it + 1
                v = v0 - (p - p_obj)/p%dv
                p = fluid%pressure(v, t)
                v0 = v
                error = abs(p - p_obj)
            end do
        else
            find_stable: block
                real(wp) :: vapor_v, liquid_v
                type(scalar_property) :: vapor_ar, liquid_ar

                vapor_v = fluid%vsolve(p_obj, t, "vapor")
                liquid_v = fluid%vsolve(p_obj, t, "liquid")

                vapor_ar = fluid%residual_helmholtz(vapor_v, t)
                liquid_ar = fluid%residual_helmholtz(liquid_v, t)

                if (vapor_ar%dt < liquid_ar%dt) then
                    v = vapor_v
                else
                    v = liquid_v
                end if
            end block find_stable
        end if
    end function solve_volume
end module fluid
