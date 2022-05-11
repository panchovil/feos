module file_operations
    !! Module to work with files
    implicit none
    integer :: out_id !! ID of an output file

    contains

    character(len=260) function outfile_with_id(str, id)
        !!From an output file name and id return a string with the style:
        !!  <name><id>.txt
        character(len=260), intent(in) :: str !! Generic name
        integer, intent(in) :: id !! Output file ID

        character(len=50) :: id_str
        write(id_str, *) id
        outfile_with_id = trim(trim(str) // adjustl(id_str)) // '.txt'
    end function
end module
