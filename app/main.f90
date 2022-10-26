program main
   use feos, only: wp, scalar_property, binary_property, PengRobinson, PR, ClassicVdW, CubicFluid
   implicit none

   integer, parameter :: n=3
   type(scalar_property) :: a(n), b(n), d, bmix, d1
   type(scalar_property) :: ar

   type(PengRobinson) :: compounds(3)
   real(wp) :: moles(n), V, T=250.0_wp

   type(ClassicVdW) :: mixing_rule
   type(binary_property) :: aij, bij
   real(wp) :: kij(n, n)
   real(wp) :: lij(n, n)

   type(CubicFluid) :: mixture

   integer :: i

   compounds(1) = PR("methane", tc=191.15_wp, pc=46.41_wp, w=0.0115_wp)
   compounds(2) = PR("ethane",  tc=305.3_wp,  pc=49.0_wp,  w=0.099_wp)
   compounds(3) = PR("propane", tc=369.9_wp,  pc=42.5_wp,  w=0.1521_wp)

   kij = reshape(&
       [0.0, 0.4, 0.3, &
        0.4, 0.0, 0.1, &
        0.3, 0.1, 0.0],&
        [n, n] &
       )
   lij = 0*kij

   mixing_rule = ClassicVdW(kij=kij, lij=lij)

   compounds%moles = 1.0_wp/n
   mixture%components = compounds
   mixture%mixing_rule = mixing_rule

   t = 250
   v = 10_wp
   ar =  mixture%residual_helmholtz(v, t)
   print *, ar%val, ar%dv, ar%dv2, ar%dtv
   print *, ar%dn
end program main

!subroutine a_b(compounds)
!   use feos
!   class(CubicEoS) :: compounds(:)
!   type(scalar_property) :: a(size(compounds)), b(size(compounds))
!
!   a = compounds%a_parameter()
!   b = compounds%b_parameter()
!   
!   print *, "a:", a%val
!   print *, "b:", b%val
!end subroutine
!
!
!subroutine mixing(compounds)
!   use feos
!
!   class(CubicEoS) :: compounds
!   real(wp) :: kij(n, n)
!   real(wp) :: lij(n, n)
!   
!   kij = reshape(&
!       [0.0, 0.4, 0.3, &
!        0.4, 0.0, 0.1, &
!        0.3, 0.1, 0.0],&
!        [n, n] &
!       )
!   lij = 0*kij
!   
!   mixing_rule = ClassicVdW(kij=kij, lij=lij)
!end subroutine

!   print *, "-----------------"
!   print *, "aij"
!   print *, "-----------------"
!   aij = aij_classic(n, a, kij)
!   do i=1,n
!      print *, aij%val(i, :)
!   end do
! 
!   print *, "-----------------"
!   print *, "bij"
!   print *, "-----------------"
!   bij = bij_classic(n, b, lij)
!   do i=1,n
!      print *, bij%val(i, :)
!   end do
!   print *, "-----------------"
!   
!   print *, "-----------------"
!   print *, "D"
!   print *, "-----------------"
!   call mixing_rule%mix(compounds, d, bmix, d1)
!   print *, "d:", d%val
!   print *, "b:", bmix%val
!   print *, "d1:", d1%val
!   print *, "-----------------"
!
!   print *, "-----------------"
!   print *, "kij"
!   print *, "-----------------"
!   do i=1,n
!      print *, kij(i, :)
!   end do
!   print *, "-----------------"
!print *, "-----------------"
!print *, "components"
!print *, "-----------------"
!do i=1,3
!   print *, compounds(i)%name
!   print *, compounds(i)%tc, compounds(i)%pc, compounds(i)%w, 0.1
!   print *, compounds(i)%ac, compounds(i)%b, compounds(i)%k
!end do
!print *, "-----------------"
