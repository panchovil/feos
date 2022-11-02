program main
   implicit none
   integer :: nin=1, nout=2
   integer :: nc, nmodel

   open(unit=nin, file='composition')
   read(1, *) nc

   call run(nc)

contains 
   subroutine run(n)
      integer :: i, j, n
      real(8) :: p, v, t, z(n), ps(500)
      
      read (1, *) (z(j), j=1, n)
      read (1, *) nmodel
      call read2PcubicNC(n, nin, nout)

      open(unit=3, file="vs")
      do i=1,500
         read(3, *) ps(i)
      end do

      z = z/sum(z) ! normalize

      t = 250.d0
      v = 0.d0

      do i = 1, 500
         p = ps(i)

         call VCALC(0, n, 2, z, T, P, V)

         print *, v, p
      end do

   end subroutine run
end program main
