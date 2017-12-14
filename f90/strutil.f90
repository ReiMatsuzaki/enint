#include "macros.fpp"
module Mod_StrUtil
  use Mod_ErrHandle
  implicit none
contains
  subroutine is_i(str, res)
    implicit none
  
    character(*) str
    double precision  :: xd
    integer :: xi  
    logical res
    integer ierr_d, ierr_i, idx
 
    double precision eps
    
    eps = 1.0D-10
    read(str, *, iostat=ierr_d) xd
    read(str, *, iostat=ierr_i) xi
    
    if(ierr_i .ne. 0) then
       ! -- failed to convert to integer --
       res = .false.
       return
    end if
    
    if(ierr_d .ne. 0) then
       ! -- something wrong --
       throw_err("somthing wrong in is_i", 1)
       return
    end if
    
    if(abs(xd-xi) > eps) then
       ! -- different value for xd and xi --
       res = .false.
       return
    end if
  
    idx = index(str, ".")
    if(idx .eq. 0) then
       res = .true.
       return
    else
       res = .false.
       return
    end if

  end subroutine is_i
  subroutine is_d(str, res)
    character(*) str
    logical res
    double precision a
    read(str, *, err=998) a
    res = .true.
    return
998 continue
    res = .false.
    return
  end subroutine is_d
  subroutine convert_i(str, a)
    character(*) str
    integer   a
    
  
    double precision d
    double precision eps
    
    eps = 1.0D-10
    read(str, *, err=999) d
    read(str, *, err=999) a
    if(abs(d-a) > eps) then
       begin_err(1)
       write(0,*) "input data may be real number"
       write(0,*) "str: ", str
       end_err();
    end if
    return
999 continue
    begin_err(1)
    write(0,*) "failed to convert to integer"
    write(0,*) "str: ", str
    end_err()
  end subroutine convert_i
  subroutine convert_d(str, a)
    character(*) str
    double precision    a

    read(str, *, err=999) a
    return
999 continue
    begin_err(1)
    write(0,*) "failed to convert to real"
    write(0,*) "str", str
    end_err()
  end subroutine convert_d
  subroutine str_split(line, sep, n, res)
    character(*), intent(in) :: line
    character(*), intent(in) :: sep
    integer, intent(out) :: n
    character(*), intent(out) :: res(:)
    integer :: i, i0

    n = 1
    i0 = 1
    do i = 1, len_trim(line)-1
       if(line(i:i).eq.trim(sep)) then
          res(n) = line(i0:i-1)
          i0 = i+1
          n = n + 1
       end if
    end do
    res(n) = line(i0:len_trim(line))
!    if(n==1) then
!       res(1) = line
!    end if
    
  end subroutine str_split
  subroutine str2vec(str, nx, xs, in_calc_xs)
    character(100), intent(in) :: str
    integer, intent(out) :: nx
    double precision, intent(out) :: xs(:)
    logical ,intent(in), optional :: in_calc_xs
    character(50) :: lines(10)
    double precision :: x0, x1
    integer :: n, i
    logical :: calc_xs

    if(present(in_calc_xs)) then
       calc_xs = in_calc_xs
    else
       calc_xs = .false.
    end if

    call str_split(str, ",", n, lines); check_err()
    if(n < 1) then
       throw_err("invalid argument", 1)
    end if
    select case(lines(1))
    case("linspace")
       if(n .ne. 4) then
          throw_err("invalid argument", 1)
       end if       
       call convert_i(lines(4), nx)
       if(get_err().ne.0) then
          begin_err(1)
          write(0,*) "invalid argument"
          write(0,*) "str:", str
          write(0,*) "line[4]:", lines(4)
          end_err()
       end if
       if(calc_xs) then
          call convert_d(lines(2), x0); check_err()
          call convert_d(lines(3), x1); check_err()
          do i = 1, nx
             xs(i) = x0 + (i-1)*(x1-x0)/(nx-1)
          end do
       end if       
    case("scalar")       
       nx = 1
       if(calc_xs) then
          call convert_d(lines(2), xs(1)); check_err()
       end if
    case default
       begin_err(1)
       write(0,*) "str: ", str
       throw_err("not supported error", 1)
    end select
  end subroutine str2vec
end module Mod_StrUtil
