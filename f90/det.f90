
! n : number of orbital
! N : number of electron
! k : orbital index
! m : number of electron up to k orbitals
module Mod_DetConst
  integer, parameter :: NANVAL=-1
end module Mod_DetConst
module Mod_DetString
  use Mod_DetConst
  implicit none
  type Obj_DetString
     integer :: num_j, N, norbs
     logical, allocatable :: o(:,:)
     integer, allocatable :: w(:,:)     ! wj(j)    : node weight
     integer, allocatable :: y(:,:)     ! yj(j)    : arc weight for occupied case
  end type Obj_DetString
contains
  subroutine DetString_new(this, N, norbs)
    type(Obj_DetString) :: this
    integer, intent(in) :: N, norbs
    integer :: k, m
    
    this%N = N
    this%norbs = norbs
    allocate(this%o(0:norbs, 0:N))
    allocate(this%w(0:norbs, 0:N))
    allocate(this%y(1:norbs, 1:N))

    this%o(:,:) = .false.
    do k = 0, this%norbs
       do m = 0, this%N
          if(correct_node(this%N,this%norbs,m,k)) then
             this%o(k,m) = .true.
          end if
       end do
    end do

    this%w(:,:) = 0
    this%y(:,:) = 0
    this%w(0,0) = 1
    do k = 1, this%norbs
       do m = 0, this%N
          if(this%o(k,m)) then
             this%w(k,m) = this%w(k-1,m)
             if(m>0) then
                this%w(k,m) = this%w(k-1,m) + this%w(k-1,m-1)
                this%y(k,m) = this%w(k-1,m)
             end if
          end if
       end do
    end do
    
  end subroutine DetString_new
  subroutine DetString_dump(this)
    type(Obj_DetString) :: this
    integer k, m

    do k = 0, this%norbs
       do m = 0, this%N
          if(this%o(k,m)) then
             write(*,'(i3, 3x)', advance='no') this%w(k,m)
          else
             write(*,'("  -", 3x)', advance='no') 
          end if          
       end do
       write(*,*) 
    end do
    
  end subroutine DetString_dump
  function DetString_index(this) result(res)
    type(Obj_DetString) :: this
    integer res
  end function DetString_index
  function correct_node(N,norbs,m,k) result(res)
    integer, intent(in) :: N,norbs,m,k
    logical :: res

    res = (m.le.k .and. (N-m).le.(norbs-k))
    
  end function correct_node
end module Mod_DetString




module Mod_DetString_old
  ! Strings which describe set of determinant in determinant based CI.
  implicit none
  type Obj_DetString
     integer :: num_j, N, norbs
     integer, allocatable :: j0s(:), j1s(:)
     integer, allocatable :: kdj(:,:)  ! kdj(j,:) : lower level connected nodes with node j
     integer, allocatable :: wj(:)     ! wj(j)    : node weight
     integer, allocatable :: yj(:)     ! yj(j)    : arc weight for occupied case
  end type Obj_DetString
contains
  subroutine DetString_new(this, N, norbs)
    type(Obj_DetString) :: this
    integer, intent(in) :: N, norbs
    
    this%N = N
    this%norbs = norbs
    allocate(this%j0s(0:norbs), this%j1s(0:norbs))

    call new_j0j1(this)
    
    this%num_j = this%j1s(norbs)
    allocate(this%kdj(this%num_j, 0:1))
    allocate(this%wj(this%num_j))
    allocate(this%yj(this%num_j))
    
  end subroutine DetString_new
  subroutine new_j0j1(this)
    type(Obj_DetString) :: this
    integer k, m
    
    this%j0s(0) = 1
    this%j1s(0) = 1
    do k = 1, this%norbs-1
       this%j0s(k) = this%j1s(k-1)+1
       this%j1s(k) = this%j1s(k-1)
       do m = 0, this%N
          if(correct_node(this%N,this%norbs,m,k)) then
             this%j1s(k) = this%j1s(k)+1
          end if
       end do
    end do
    this%j0s(this%norbs) = this%j1s(this%norbs-1)+1
    this%j1s(this%norbs) = this%j0s(this%norbs)

  end subroutine new_j0j1
  function correct_node(N,norbs,m,k) result(res)
    integer, intent(in) :: N,norbs,m,k
    logical :: res

    res = (m.le.k .and. (N-m).le.(norbs-k))
    
  end function correct_node
end module Mod_DetString_Old
