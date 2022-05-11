module datatypes
    use constants
    implicit none

    type mix
        integer :: nmodel  !! Model to use
        integer :: nc      !! Number of components
        integer :: ntdep   !! Temperature dependence 
        integer :: ncomb   !! Combining rule

        character(len=:), dimension(:), allocatable :: names !! Components names

        real(wp), dimension(:), allocatable :: n !! Number of moles

        ! Properties
        real(wp), dimension(:), allocatable :: tc !! Critical Temperatures
        real(wp), dimension(:), allocatable :: pc !! Critical Pressures 
        real(wp), dimension(:), allocatable :: dc !! Critical Densities (from EOS)
        real(wp), dimension(:), allocatable :: w !! Accentric factors
        ! EOS parameters
        real(wp), dimension(:), allocatable :: ac !! EOS atractive parameter
        real(wp), dimension(:), allocatable :: b !! EOS repulsive parameter
        real(wp), dimension(:), allocatable :: del1 !! EOS delta_1
        real(wp), dimension(:), allocatable :: k !! k parameter to calculate the a parameter

        ! Mixing parameters
        real(wp), dimension(:, :), allocatable :: kij !! Kij matrix
        real(wp), dimension(:, :), allocatable :: kij0 !! Kij standard
        real(wp), dimension(:, :), allocatable :: kij_inf !! Kij at infinite temperature
        real(wp), dimension(:, :), allocatable :: T_star !! Reference temperature for temperature dependent Kij

        real(wp), dimension(:, :), allocatable :: lij !! lij matrix
        real(wp), dimension(:, :), allocatable :: aij !! EOS atractive parameter matrix
    end type mix

    type compound
        !! 
        real(wp) :: tc
        real(wp) :: pc
        real(wp) :: ac
        real(wp) :: b
        real(wp) :: w
        real(wp) :: k
    end type compound
end module datatypes

module mixture
    !! Module to represent a mixture of fluids, it's used to save the 
    !! multiple component's properties to be used by different subroutines.
    !! It also includes a set of subroutines to read the data and save it
    !! in the module.

    !! System of units: This units are asumed, if the user wants to use another
    !! system, the RGAS constant should be changed at the module `constants`
    !!
    !! - Volume: Liter
    !! - Pressure: bar
    !! - Temperature: Kelvin
    use constants
    implicit none

    integer :: nmodel  !! Model to use
    integer :: nc      !! Number of components
    integer :: ntdep   !! Temperature dependence 
    integer :: ncomb   !! Combining rule

    character(len=:), dimension(:), allocatable :: names !! Components names

    real(wp), dimension(:), allocatable :: n !! Number of moles

    ! Properties
    real(wp), dimension(:), allocatable :: tc !! Critical Temperatures
    real(wp), dimension(:), allocatable :: pc !! Critical Pressures 
    real(wp), dimension(:), allocatable :: dc !! Critical Densities (from EOS)
    real(wp), dimension(:), allocatable :: w !! Accentric factors
    ! EOS parameters
    real(wp), dimension(:), allocatable :: ac !! EOS atractive parameter
    real(wp), dimension(:), allocatable :: b !! EOS repulsive parameter
    real(wp), dimension(:), allocatable :: del1 !! EOS delta_1
    real(wp), dimension(:), allocatable :: k !! k parameter to calculate the a parameter

    ! Mixing parameters
    real(wp), dimension(:, :), allocatable :: kij !! Kij matrix
    real(wp), dimension(:, :), allocatable :: kij_0 !! Kij standard
    real(wp), dimension(:, :), allocatable :: kij_inf !! Kij at infinite temperature
    real(wp), dimension(:, :), allocatable :: T_star !! Reference temperature for temperature dependent Kij

    real(wp), dimension(:, :), allocatable :: lij !! lij matrix
    real(wp), dimension(:, :), allocatable :: aij !! EOS atractive parameter matrix

    contains

    subroutine setup(nin, filename)
        !! This subroutine will be used to read data files and get the
        !! mixture properties
        integer, intent(in) :: nin
        character(len=:), allocatable, intent(in) :: filename

        read(nin, *) nc
        read(nin, *) nmodel
        read(nin, *) ncomb, ntdep
    end subroutine setup
end module mixture
