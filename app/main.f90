program main
   use constants
   use properties
   use cubic_eos
   use mixing_rules

   implicit none

   integer, parameter :: n=3
   type(PengRobinson) :: compounds(n)

   type(scalar_property) :: a(n), d, b, d1
   type(ClassicVdW) :: mixing_rule

   real(wp) :: kij(n, n)
   real(wp) :: lij(n, n)

   integer :: i

   compounds(1) = PR("methane", tc=191.15_wp, pc=46.41_wp, w=0.0115_wp)
   compounds(2) = PR("ethane",  tc=305.3_wp,  pc=49.0_wp,  w=0.099_wp)
   compounds(3) = PR("propane", tc=369.9_wp,  pc=42.5_wp,  w=0.1521_wp)

   compounds(1)%moles = 1.0_wp/3.0_wp
   compounds(2)%moles = 1.0_wp/3.0_wp
   compounds(3)%moles = 1.0_wp/3.0_wp

   kij = reshape(&
       [0.0, 0.4, 0.3, &
        0.4, 0.0, 0.1, &
        0.3, 0.1, 0.0],&
        [n, n] &
       )
   lij = 0*kij

   mixing_rule = ClassicVdW(kij=kij, lij=lij)
   call mixing_rule%mix(compounds, d, b, d1)

   print *, d%val, b%val, d1%val
end program main
