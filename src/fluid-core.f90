module fluid_core
    !! Abstract types related to fluid derived types.
    type, abstract :: BaseFluid
        !! Abstract Fluid
    contains
        !> Function that calculates the residual Helmholtz Energy
        procedure(vt_property), deferred :: residual_helmholtz
    end type BaseFluid
    
    abstract interface 
    elemental function vt_property(self, v, t) result(prop)
        !! Interface of a type procedure that receives volume and temperature
        !! and return some kind of scalar property.
        use properties, only: scalar_property
        use constants
        import BaseFluid
        class(BaseFluid), intent(in) :: self

        real(wp), intent(in) :: v
        real(wp), intent(in) :: t
        type(scalar_property) :: prop

    end function
    end interface
end module fluid_core
