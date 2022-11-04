program main
   use feos, only: wp, scalar_property, binary_property, PengRobinson, PR, ClassicVdW, CubicFluid, R
   use root_module
   implicit none

   integer, parameter :: n=3
   real(wp), parameter:: RGAS=0.08314472
   real(wp) :: preal, ps(1000)
   type(scalar_property) :: a(n), b(n), d, bmix, d1
   type(scalar_property) :: ar, p, lnfug(n)

   type(PengRobinson) :: compounds(n)
   real(wp) :: moles(n), v, vsolved, t

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

   t = 150_wp
   v = 5_wp
   ar = mixture%residual_helmholtz(v, t)   
   print *, ar%val

   do i = 150, 800
      ps(i) = real(i, 8)/1000

      p%val = ps(i)

      vsolved = mixture%vsolve(p%val, t, max_it=50)
      ar = mixture%residual_helmholtz(vsolved, t)

      print *, vsolved, p%val, ar%val
   end do

end program main
