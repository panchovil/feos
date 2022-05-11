program main
   use json_module
   implicit none

   type(json_file) :: json_data

   logical :: found
   integer :: n, i

   real(8), allocatable :: z(:)
   real(8), allocatable :: kijs(:, :)
   real(8), allocatable :: kij(:)

   character*50 :: id

   character(len=:), allocatable :: name

   call json_data%initialize()
   call json_data%load("mixfile.json")
   call json_data%print()

   call json_data%get('z', z, found)
   n = size(z)

   do i = 1, n
      write (id, *) i
      id = trim('compounds('//trim(adjustl(id))//').name')
      print *, id
      call json_data%get(id, name, found)
      print *, name
   end do
end program main
