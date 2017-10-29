! #include "macros.fpp"
!     Error handle utilities
module mod_ErrHandle
  implicit none  
  logical show_message_q
  integer ierr
contains
  subroutine ErrHandle_new
    show_message_q = .true.
    ierr = 0
  end subroutine ErrHandle_new
  subroutine ErrHandle_delete
    ierr = 0
  end subroutine ErrHandle_delete
  subroutine ErrHandle_err(msg, ierr0, file, line)
    character(*) msg
    integer ierr0
    character(*) file
    integer line
    if(show_message_q) then
       write(0, '(A, ":", I0, ": ", A)') file, line, msg
    end if
    ierr = ierr0
  end subroutine ErrHandle_err
  subroutine open_w(ifile, filename)
    integer, intent(in) ::  ifile
    character(*), intent(in) ::  filename
    
    open(ifile, file=filename, status="replace", err=999)
    return
    
999 continue
    call err_with_file_line("failed to open file", 1, __FILE__, __LINE__)
    call err_ss("filename: ", filename)
    return
  end subroutine open_w
  subroutine open_r(ifile, filename)
    integer ifile
    character(*) filename
    
    open(ifile, file=filename, status='old',err=999)
    return
    
999 continue
    call err_with_file_line("failed to open file", 1, __FILE__, __LINE__)
    call err_ss("filename: ", filename)
    return
    
  end subroutine open_r
end module mod_ErrHandle

