#include "macros.fpp"

module Mod_ArgParser
  use Mod_ErrHandle  
  use Mod_fjson
  implicit none
  integer :: num_
  character(100), allocatable :: names_(:)
  integer, allocatable :: itypes_(:)
  integer, allocatable :: nums_(:)
  integer :: idx_
  type(object) :: args_  
contains
  subroutine ArgParser_new
    idx_ = 1
    num_ = 10
    call alloc(num_)
  end subroutine ArgParser_new
  subroutine ArgParser_delete
    deallocate(names_)
    deallocate(itypes_)
    deallocate(nums_)
  end subroutine ArgParser_delete
  subroutine ArgParser_add(name, itype, num)
    character(*), intent(in) :: name
    integer, intent(in) :: itype, num

    names_(idx_) = name
    itypes_(idx_) = itype
    nums_(idx_) = num
    idx_ = idx_ + 1
    
  end subroutine ArgParser_add
  subroutine ArgParser_parse(o)
    use Mod_math, only : is_i, convert_i, convert_d
    type(object) :: o
    integer :: narg
    integer :: iarg, idx, num, itype, j, vali
    double precision :: vald
    character :: ele*100
    character(100) :: opts_name
    logical :: numq

    narg = iargc()
    iarg = 0
    do 
       iarg = iarg+1
       if(iarg > narg) then
          exit
       end if

       call getarg(iarg, ele)
       
       if(ele(1:1).eq."-") then
          call is_i(ele(2:2), numq)
          if(numq) then
             throw_err("-(numeric) is not supported", 1)
          end if
          if(ele(2:2).eq."-") then
             throw_err("-- type option is not supported", 1)
          end if
          opts_name = ele(2:100)
          idx = get_idx(trim(opts_name))

          itype = itypes_(idx)
          num = nums_(idx)
          if(num.ne.1) then
             throw_err("only num==1 is supported", 1)
          end if
          do j = 1, num
             iarg = iarg+1
             call getarg(iarg, ele)
             select case(itype)
             case(TYPE_I)
                call convert_i(ele, vali); check_err()
                call object_set_i(o, opts_name, vali)
             case(TYPE_S)
                call object_set_s(o, opts_name, ele)
             case(TYPE_D)
                call convert_d(ele, vald); check_err()
                call object_set_d(o, opts_name, vald)
             case default
                throw_err("unsupported TYPE", 1)
             end select
          end do
       end if
    end do
    
  end subroutine ArgParser_parse
  function get_idx(name) result(res)
    character(*), intent(in) :: name
    integer res
    integer i
    res = 0
    do i = 1, num_
       if(names_(i).eq.name) then
          res = i
          return
       end if
    end do
    throw_err("failed to find name", 1)
  end function get_idx
  subroutine alloc(num)
    integer, intent(in) :: num
    if(allocated(nums_)) then
       if(size(nums_) < num) then
          deallocate(names_)
          deallocate(itypes_)
          deallocate(nums_)
       end if
    end if

    if(.not. allocated(nums_)) then
       allocate(names_(num))
       allocate(itypes_(num))
       allocate(nums_(num))
    end if
          
  end subroutine alloc
end module Mod_ArgParser
