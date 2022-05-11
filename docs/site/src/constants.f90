module constants
   implicit none

   integer, parameter :: wp = 8
   real(wp), parameter :: RGAS = 0.08314472d0
   real(wp), parameter :: A0 = 0.0017, B0 = 1.9681, C0 = -2.7238
   real(wp), parameter :: A1 = -2.4407, B1 = 7.4513, C1 = 12.504
   real(wp), dimension(6) :: D = [0.428363, 18.496215, 0.338426, &
                                  0.660, 789.723105, 2.512392]
   real(wp), parameter :: ERRMAX = 1.D-8
end module constants
