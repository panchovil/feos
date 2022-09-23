module properties
   !! Datatype to represent different kind of properties with the objective
   !! to ease accesibility to derivatives and don't use an excess of arguments.
   use constants, only: wp, R

   implicit none

   type :: scalar_property(n)
       !! Scalar propery
       integer, len :: n
       real(wp) :: val !! Value

       real(wp) :: dt        !! First deritvative with temperature
       real(wp) :: dt2       !! Second derivative with temperature
       real(wp) :: dtn(n)    !! Cross derivative with composition and temperature

       real(wp) :: dv        !! First deritvative with volume
       real(wp) :: dv2       !! Second derivative with volume
       real(wp) :: dvn(n)    !! Cross derivative with composition and volume

       real(wp) :: dtv       !! Cross derivative with temperature and volume

       real(wp) :: dn(n)     !! First derivative with composition
       real(wp) :: dn2(n, n) !! Second derivative with composition
   end type scalar_property

   type :: binary_property(n)
      !! Binary property
      integer, len :: n
      real(wp) :: val(n, n)  !! Value

      real(wp) :: dt(n, n)   !! First deritvative with temperature
      real(wp) :: dt2(n, n)  !! Second derivative with temperature

      real(wp) :: dv(n, n)   !! First deritvative with volume
      real(wp) :: dv2(n, n)  !! Second derivative with volume

      real(wp) :: dtv(n, n) !! Crossed derivative with temperature and volume
   end type
end module properties

