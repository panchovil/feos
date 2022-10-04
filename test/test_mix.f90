program main
    !use cubic_eos
    !use mixing_rules
    !use multicomponent_eos

    implicit none

    !type(srk) :: methane, nitrogen, ethane
    !type(srk) :: compounds(3)
    !type(kij_exp_t) :: mixing_rule
    !type(mixture) :: fluid

    !real(8), dimension(3, 3) :: kij, kij_0, kij_inf, T_star, dkijdt, dkij2dt2, lij

    !character(len=*), parameter :: fmt1 = '(F15.4, F15.4, F15.4)'

    !integer :: i

    !! -------------------------------------------------------------------------
    !! Reading "input"
    !methane = srk("methane", 0.01, 2, 3, 4, 5, 6)
    !nitrogen = srk("nitrogen", 0.07, 8, 9, 10, 11, 12)
    !ethane = srk("ethane", 0.019, 5, 2, 15, 1, 2)

    !kij = reshape([&
    !    0.0, 0.15, 0.20, &
    !    0.15, 0.0, 0.25, &
    !    0.20, 0.25, 0.0  &
    !], [3, 3])
    !dkijdt = 0*kij
    !dkij2dt2 = 0*kij

    !kij_0 = kij/5.d0
    !kij_inf = kij/10.d0
    !T_star = kij/kij * 300
    !
    !lij = reshape([0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 3])

    !! -------------------------------------------------------------------------
    !! Defining objets
    !!mixing_rule%kij = kij
    !mixing_rule%lij = lij

    !mixing_rule%kij_0   = kij_0
    !mixing_rule%kij_inf = kij_inf
    !mixing_rule%T_star  = T_star

    !compounds = [methane, nitrogen, ethane]

    !fluid% T = 250.d0
    !fluid% compounds = compounds
    !fluid% mixing_rule = mixing_rule
    !fluid% concentrations = [1.d0, 2.d0, 3.d0]

    !! -------------------------------------------------------------------------
    !! Running methods
    !call fluid%mix()
    !! -------------------------------------------------------------------------

    !print *, "-----------------------------------------------------------------"
    !print *, (trim(fluid%compounds(i)%name), i=1,size(compounds))
    !print *, "a parameters-----------------------------------------------------"
    !print fmt1, (fluid%compounds(i)%a, i=1,size(compounds))
    !print *, "b parameters-----------------------------------------------------"
    !print fmt1, (fluid%compounds(i)%b, i=1,size(compounds))

    !associate(mr => fluid%mixing_rule)
    !
    !select type(mr)
    !class is (quadratic)
    !print *, "kij matrix-------------------------------------------------------"
    !do i = 1, size(compounds)
    !    print fmt1, mr%kij(i, :size(compounds))
    !end do

    !print *, "aij matrix-------------------------------------------------------"
    !do i = 1, size(compounds)
    !    print fmt1, mr%aij(i, :size(compounds))
    !end do

    !print *, "bij matrix-------------------------------------------------------"
    !do i = 1, size(compounds)
    !    print fmt1, mr%bij(i, :size(compounds))
    !end do
    !end select

    !end associate
end program main
