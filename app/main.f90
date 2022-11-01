program main
   use feos, only: wp, scalar_property, binary_property, PengRobinson, PR, ClassicVdW, CubicFluid, R
   use root_module
   implicit none

   integer, parameter :: n=3
   real(8), parameter:: RGAS=0.08314472
   real(8) :: preal
   type(scalar_property) :: a(n), b(n), d, bmix, d1
   type(scalar_property) :: ar, p, lnfug(n)

   type(PengRobinson) :: compounds(n)
   real(wp) :: moles(n), v, vsolved, t=250.0_wp

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

   ar = mixture%residual_helmholtz(v, t)
   print *, ar%val

   do i = 1, 500
      v = real(i, 8)/200
      ar = mixture%residual_helmholtz(v, t)
      preal = RGAS * t/v - ar%dv
      vsolved = mixture%vsolve(preal, t, max_it=25)

      print *, preal, v, vsolved, v - vsolved
   end do
end program main
