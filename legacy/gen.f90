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
      real(8) :: t=250.d0, aij(n,n), aijdt(n,n), daijdt2(n,n), daijdt(n,n), z(n)
      real(8) :: d, ddi(n), ddit(n), ddij(n,n), dddt, dddt2
      real(8) :: Ar, Arv, artv, arv2, arn(n), arvn(n), artn(n), arn2(n, n)

      real(8) :: temps(3) = [200, 250, 300]
      
      D = 0.d0
      ddi = 0
      ddit = 0
      ddij = 0
      dddt = 0
      dddt2 = 0

      read (1, *) (z(j), j=1, n)
      read (1, *) nmodel
      call read2PcubicNC(n, nin, nout)

      call aijTder(1, n, T, aij, daijdt, daijdt2)
      print *, "--------------"
      print *, "aij"
      do i=1,n
         print *, aij(i, :)
      end do
      print *, "--------------"
      print *, "aij dt"
      do i=1,n
         print *, daijdt(i, :)
      end do
      print *, "--------------"
      print *, "aij dt2"
      do i=1,n
         print *, daijdt2(i, :)
      end do
      print *, "--------------"


      call DandTnder(1, n, T, z, D, dDi, dDiT, dDij, dDdT, dDdT2)

      print *, "D"
      print *, D, ddi(:), shape(ddi)
      print *, "--------------"

      z = z/sum(z) ! normalize

      do i = 1, 3
         T = 250 ! temps(i)
         call HelmSRKPR(n, 2, 1, z, 10.0, T, Ar, ArV, ArtV, Arv2, Arn, arVn, ArTn, Arn2)
         print *, Ar, arv, arv2, artv
         print *, arn

         do j=1,3
            print *, arn2(:, j)
         end do
         
         stop
      end do
   end subroutine run
end program main
