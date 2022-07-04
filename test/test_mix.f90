program main
    use cubic_eos
    use multicomponent_eos

    implicit none

    type(srk) :: methane, nitrogen, ethane
    type(srk) :: compounds(3)
    type(mixture) :: fluid

    real(8), dimension(3, 3) :: kij, kij_0, kij_inf, T_star, dkijdt, dkij2dt2, lij

    class(kij_exp_t), allocatable :: mixing_rule
    character(len=*), parameter :: fmt1 = '(F15.4, F15.4, F15.4)'

    integer :: i

    ! -------------------------------------------------------------------------
    ! Reading "input"
    methane = srk("methane", 1, 2, 3, 4, 5, 6)
    nitrogen = srk("nitrogen", 7, 8, 9, 10, 11, 12)
    ethane = srk("ethane", 9, 5, 2, 15, 1, 2)

    kij = reshape([0 , 15, 20, 15, 0, 25, 20, 25, 0], [3, 3])
    dkijdt = 0*kij
    dkij2dt2 = 0*kij

    kij_0 = kij/5.d0
    kij_inf = kij/10.d0
    T_star = kij*25.d0
    
    lij = reshape([0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 3])

    ! -------------------------------------------------------------------------
    ! Defining objets
    mixing_rule = kij_exp_t(kij, lij, lij, lij, dkijdt, dkij2dt2, kij_0, kij_inf, T_star)
    compounds = [methane, nitrogen, ethane]

    fluid%compounds = compounds
    fluid%mixing_rule = mixing_rule

    ! -------------------------------------------------------------------------
    ! Running methods
    call fluid%mixing_rule%mix(250.d0, fluid%compounds)
    ! -------------------------------------------------------------------------

    print *, "a parameters-----------------------------------------------------"
    print *, (trim(fluid%compounds(i)%name), i=1,size(compounds))
    print *, "a parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%a, i=1,size(compounds))
    print *, "b parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%b, i=1,size(compounds))

    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%kij(i, :size(compounds))
    end do

    print *, "aij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%aij(i, :size(compounds))
    end do

    print *, "bij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%bij(i, :size(compounds))
    end do


    !call fluid%mixing_rule%mix(500.d0, fluid%compounds)

    print *, "a parameters-----------------------------------------------------"
    print *, (trim(fluid%compounds(i)%name), i=1,size(compounds))
    print *, "a parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%a, i=1,size(compounds))
    print *, "b parameters-----------------------------------------------------"
    print fmt1, (fluid%compounds(i)%b, i=1,size(compounds))

    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%kij(i, :size(compounds))
    end do

    print *, "aij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%aij(i, :size(compounds))
    end do

    print *, "bij matrix-------------------------------------------------------"
    do i = 1, size(compounds)
        print fmt1, fluid%mixing_rule%bij(i, :size(compounds))
    end do
end program main
