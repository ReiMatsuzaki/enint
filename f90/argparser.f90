module Mod_ArgParser
  use Mod_ErrHandle  
  use Mod_fjson
  implicit none
  character(100), allocatable :: short_names_(:)
  character(100), allocatable :: long_names_(:)
  integer, allocatable :: itypes_(:)
  integer, allocatable :: nums_(:)
  integer :: idx_
  type(object) :: args_  
contains
  subroutine ArgParser_new
    idx_ = 1
    call alloc(10)
  end subroutine ArgParser_new
  subroutine ArgParser_delete
    deallocate(short_names_)
    deallocate(long_names_)
    deallocate(itypes_)
    deallocate(nums_)
  end subroutine ArgParser_delete
  subroutine ArgParser_add(short_name, long_name, itype, num)
    character(*), intent(in) :: short_name, long_name
    integer, intent(in) :: itype, num

    short_names_(idx_) = short_name
    long_names_(idx_) = long_name
    itypes_(idx_) = itype
    nums_(idx_) = num
    idx_ = idx_ + 1
    
  end subroutine ArgParser_add
  subroutine ArgParser_parse(o)
    type(object) :: o
    integer :: narg
    integer :: i
    character :: ele*100

    narg = iargc()
    do i = 1, narg
       call getarg(i, ele)
       write(*,'(a)') ele(1:2)
    end do
    
  end subroutine ArgParser_parse
  subroutine alloc(num)
    integer, intent(in) :: num
    if(allocated(nums_)) then
       if(size(nums_) < num) then
          deallocate(short_names_)
          deallocate(long_names_)
          deallocate(itypes_)
          deallocate(nums_)
       end if
    end if

    if(.not. allocated(nums_)) then
       allocate(short_names_(num))
       allocate(long_names_(num))
       allocate(itypes_(num))
       allocate(nums_(num))
    end if
          
  end subroutine alloc
end module Mod_ArgParser
