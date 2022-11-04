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
      real(8) :: p, v, t, z(n)
      real(8) :: Ar, ArV, ArTV, ArV2, Arn(n), ArVn(n), ArTn(n), Arn2(n, n)
      real(8), parameter :: R=0.08314472
      
      read (1, *) (z(j), j=1, n)
      read (1, *) nmodel
      call read2PcubicNC(n, nin, nout)

      z = z/sum(z) ! normalize

      t = 150.d0
      v = 5.d0
      
!      v = 25.05227907637569
!      call ArVnder(n, 2, 1, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
!      print *, ar, arv
!
!      v = 0.050351279383939546
!      call ArVnder(n, 2, 1, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
!      print *, ar, arv
!
!      call exit

      do i = 150, 800
         p = real(i, 8)/1000
         call VCALC(0, n, 2, z, T, P, V)
         call ArVnder(n, 2, 1, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
         if (i .eq. 488 .or. i .eq. 489) then
            print  *, (Ar + V*P)/(T*R) - LOG(V)
         end if
         print *, v, p, ar 
      end do

      p = sum(z) * R * t/v - arv
      call VCALC(0, n, 2, z, T, P, V)
   end subroutine run
end program main
