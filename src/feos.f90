module feos
    use constants
    use properties
    use cubic_eos
    use mixing_rules
    use fluid

    implicit none
    private

    ! Cubic Equations of state
    public :: PengRobinson, PR
    public :: SoaveRedlichKwong, SRK
    public :: BaseCEOS, CubicEoS
    ! Properties
    public :: scalar_property, binary_property
    ! Constnats
    public :: wp, R
    ! Mixing rules
    public :: ClassicVdW
    public :: CubicFluid
end module feos
