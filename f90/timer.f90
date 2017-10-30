#include "macros.fpp"

module Mod_Timer
  use Mod_ErrHandle
  implicit none
  private
  type Record
     character(100) :: key
     double precision :: tsum
     integer :: t0
     integer :: num
     logical :: act
  end type Record
  character(10) :: name_
  integer :: size_, capacity_
  type(Record), allocatable :: rec_(:)
  integer t0_
  logical print_key_
  public Timer_new, Timer_delete, Timer_time, Timer_begin, Timer_end, Timer_result
contains
  subroutine init(num)
    integer, intent(in) :: num
    integer i
    allocate(rec_(num))
    do i = 1, num
       rec_(i)%key = ""
       rec_(i)%tsum = 0
       rec_(i)%t0 = 0
       rec_(i)%num = 0
       rec_(i)%act = .false.
    end do
  end subroutine init
  subroutine extend
    integer num0, num1
    type(Record), allocatable :: tmp(:)
    
    num0 = capacity_
    num1 = 2*num0

    allocate(tmp(num0))
    tmp(:) = rec_(:)

    deallocate(rec_)
    call init(num1); check_err()
    rec_(:num0) = tmp(:)
    
    deallocate(tmp)
    capacity_ = num1
    
  end subroutine extend
  function find(key) result(res)
    character(*), intent(in) :: key
    integer :: res
    integer i

    res = 0
    do i = 1, size_
       if(trim(rec_(i)%key)==trim(key)) then
          res = i
          return
       end if
    end do
    
  end function find
  function Timer_time(key) result(res)
    character(*), intent(in) :: key
    double precision res
    integer idx

    idx = find(trim(key))

    res = rec_(idx)%tsum
    
  end function Timer_time
  subroutine Timer_new(name, print_key)
    character(*), intent(in) :: name
    logical, intent(in) :: print_key
    integer, parameter :: NUM0=1024
    name_ = name
    size_ = 0
    capacity_ = NUM0
    print_key_ = print_key
    call system_clock(t0_)
    call init(NUM0)
  end subroutine Timer_new
  subroutine Timer_delete
    deallocate(rec_)
  end subroutine Timer_delete
  subroutine Timer_begin(key)
    character(*), intent(in) :: key
    integer idx
    if(print_key_) then
       write(*,'(A,":",A," begin")') trim(name_), key
    end if
    idx = find(key)
    if(idx==0) then
       if(size_+1>capacity_) then
          call extend
       end if       
       size_ = size_+1
       rec_(size_)%key = key
       rec_(size_)%tsum = 0
       call system_clock(rec_(size_)%t0)
       rec_(size_)%num = 0
       rec_(size_)%act = .true.
    else
       if(rec_(idx)%act) then
          begin_err(1)
          write(0,*) "Timer_start is called for active timer key."
          write(0,*) "key:", key
          end_err()
       end if
       call system_clock(rec_(idx)%t0)
       rec_(idx)%act = .true.
    end if
  end subroutine Timer_begin
  subroutine Timer_end(key)
    character(*), intent(in) :: key
    integer idx
    integer dt, t0, t1, t_rate, t_max

    idx = find(key)
    if(idx==0) then
       begin_err(1)
       write(0,*) "Timer_end is called for not started record"
       write(0,*) "key:", key
       end_err()
    else
       if(.not.rec_(idx)%act) then
          begin_err(1)
          write(0,*) "Timer_end is called for non active started record"
          write(0,*) "key:", key
          end_err()
       end if
       t0 = rec_(idx)%t0
       call system_clock(t1, t_rate, t_max); check_err()
       if(t1<t0) then
          dt = (t_max-t0) + t1 + 1
       else
          dt = t1-t0
       end if
       rec_(idx)%tsum = rec_(idx)%tsum + (1.0d0*dt)/t_rate
       rec_(idx)%num = rec_(idx)%num+1
       rec_(idx)%act = .false.
    end if

    if(print_key_) then
       write(*,'(A,":",A," end")') trim(name_), key
    end if
    
  end subroutine Timer_end
  subroutine Timer_result(in_ifile)
    integer, optional :: in_ifile
    integer i, ifile
    integer t1, t_rate, t_max, dt

    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if
    
    write(ifile,'(A,":Time_result begin")') trim(name_)

    call system_clock(t1, t_rate, t_max); check_err()
    if(t1<t0_) then
       dt = (t_max-t0_) + t1 + 1
    else
       dt = t1-t0_
    end if
    write(ifile,*) "total time;", (1.0d0*dt)/t_rate
       
    do i = 1, size_
       if(rec_(i)%act) then
          begin_err(1)
          write(0,*) "active record is found."
          write(0,*) "key:", rec_(i)%key
          end_err()
       end if
       write(ifile,'(A30, I10, f20.5)') trim(rec_(i)%key), rec_(i)%num, rec_(i)%tsum
    end do
    write(ifile,'(A,":Time_result end")') trim(name_)
  end subroutine Timer_result
end module Mod_timer
