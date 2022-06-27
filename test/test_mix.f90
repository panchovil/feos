program main
    use cubic_eos
    use multicomponent_eos

    type(pure_compound) :: methane, nitrogen, ethane
    type(pure_compound) :: compounds(3)
    type(mixture) :: fluid
    real(8), dimension(3, 3) :: kij, kij_0, kij_inf, T_star, dkijdt, dkij2dt2

    class(kij_constant), allocatable :: kij_calculator
    type(lij_constant), allocatable :: lij_calculator
    character(len=*), parameter :: fmt1 = '(F15.4, F15.4, F15.4)'

    methane = pure_compound("methane", 1, 2, 3, 4, 5, 6)
    nitrogen = pure_compound("nitrogen", 7, 8, 9, 10, 11, 12)
    ethane = pure_compound("ethane", 9, 5, 2, 15, 1, 2)

    kij = reshape([0 , 15, 20, 15, 0, 25, 20, 25, 0], [3, 3])
    dkijdt = 0*kij
    dkij2dt2 = 0*kij

    kij_0 = kij/5.d0
    kij_inf = kij/10.d0
    T_star = kij*25.d0
        
    kij_calculator = kij_exp_t(kij, dkijdt, dkij2dt2, kij_0, kij_inf, T_star)
    lij_calculator = lij_constant(reshape([0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 3]))

    compounds = [methane, nitrogen, ethane]
    fluid = mixture(compounds, kij_calculator, lij_calculator)

    call fluid%mix(250.d0)

    print *, "a parameters-----------------------------------------------------"
    print *, (trim(fluid%compounds(i)%name), i=1,size(compounds))
    print *, "a parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%a, i=1,size(compounds))
    print *, "b parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%b, i=1,size(compounds))

    do i = 1, size(compounds)
        print fmt1, kij_calculator%kij(i, :size(compounds))
    end do

    print *, "aij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%aij(i, :)
    end do

    print *, "bij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%bij(i, :)
    end do
end program main
