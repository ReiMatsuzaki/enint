#include "macros.fpp"

! i : level
! j : DR number

module Mod_DRTConst
  integer, parameter :: NANVAL=-1
  integer :: delta_abc_(0:3,3) = reshape((/0,0,1,1,  0, 1,-1, 0,  1, 0, 1, 0/), (/4,3/))
end module Mod_DRTConst

module Mod_DRTBuilder
  ! builder for Distinct Row Table(DRT).

  use Mod_DRTConst, only : NANVAL
  implicit none
  integer :: nj_ ! number of j
  integer, allocatable :: j0s_(:), j1s_(:) ! j0s_(i)/j1s_(i) is min/max of j.
  integer, allocatable :: j_table_(:,:,:)  ! j_table_(a,b,c) return DR number j.
  integer, allocatable :: abc_table_(:,:)  ! abc_table(j,:) return (a,b,c) for j.
  integer :: current_j_
contains
  subroutine DRTBuilder_new(N, TwoS, norbs)
    !
    ! Inputs
    ! -------
    ! N    : number of electron
    ! TwoS : twice of total spin quantum number
    ! norbs: number of orbital
    
    integer, intent(in) :: N, TwoS, norbs
    integer :: na, nb, nc, abc(3), i, j0, j1
    
    ! -- compute table --
    current_j_ = 1
    na = int(N/2)
    nb = int(TwoS + norbs)
    nc = int(norbs)

    allocate(j0s_(0:norbs), j1s_(0:norbs))
    allocate(j_table_(0:na, 0:nb, 0:nc))
    allocate(abc_table_((na+1)*(nb+1)*(nc+1),3))
    j_table_(:,:,:) = NANVAL
    abc_table_(:,:) = NANVAL

    abc(1) = (N-TwoS)/2
    abc(2) = TwoS
    abc(3) = norbs-abc(1)-abc(2)
    call get_j(abc, j0s_(norbs))
    j1s_(norbs) = j0s_(norbs)
    do i = norbs, 1, -1
       call next_j0j1(j0s_(i), j1s_(i), j0, j1)
       j0s_(i-1) = j0
       j1s_(i-1) = j1
    end do

    nj_ = j1s_(0)
    
  end subroutine DRTBuilder_new
  subroutine DRTBuilder_delete
    deallocate(j0s_, j1s_)
    deallocate(j_table_, abc_table_)
  end subroutine DRTBuilder_delete
  subroutine get_j(abc, j)
    ! return DR number j for abc with registration of (j,abc) pair
    ! Inputs
    ! ------
    ! abc : (a,b,c)
    !
    ! Returns
    ! -------
    ! j : DR number
    
    integer, intent(in) :: abc(3)
    integer, intent(out) :: j
    integer :: a, b, c
    a = abc(1)
    b = abc(2)
    c = abc(3)
    if(j_table_(a,b,c) .eq. NANVAL) then
       j = current_j_
       j_table_(a,b,c) = j
       abc_table_(j,:) = abc(:)
       current_j_ = j + 1
    else
       j = j_table_(a,b,c)
    end if
  end subroutine get_j
  subroutine next_j0j1(j0,j1, jj0, jj1)
    ! from (j0,j1) of i level, compute (jj0,jj1) of i-1 level.
    
    use Mod_DRTConst, only : delta_abc_
    integer, intent(in) :: j0, j1
    integer, intent(out) :: jj0, jj1
    integer :: j, abc(3), abc2(3), d, jj

    jj0 = j1+1
    do j = j0, j1
       abc = abc_table_(j,:)
       do d = 0, 3
          abc2 = abc - delta_abc_(d,:)
          if(correct_abc(abc2)) then
             call get_j(abc2, jj)             
          end if
       end do
    end do
    jj1 = jj
    
  end subroutine next_j0j1
  subroutine DRTBuilder_abck(abcj, kdj, ukdj)
    use Mod_DRTConst, only : delta_abc_
    integer, intent(out) :: abcj(:,:), kdj(:,0:), ukdj(:,0:)
    integer :: j, d, abc_d(3), jj

    kdj(:,:) = NANVAL
    ukdj(:,:) = NANVAL    
    do j = 1, nj_
       abcj(j,:) = abc_table_(j,:)
       do d = 0, 3
          abc_d(:) = abcj(j,:) - delta_abc_(d,:)
          if(correct_abc(abc_d)) then
             call get_j(abc_d, jj)
             kdj(j,d) = jj
             ukdj(jj,d) = j
          end if
       end do
    end do
    
  end subroutine DRTBuilder_abck
  subroutine DRTBuilder_xy(kdj, ukdj, xj, ydj, uxj, uydj)
    integer, intent(in) :: kdj(:,0:), ukdj(:,0:)
    integer, intent(out) :: xj(:), ydj(:,0:), uxj(:), uydj(:,0:)
    integer j, dd, d
    
    xj(nj_) = 1
    do j = nj_, 1, -1
       dd = NANVAL
       do d = 0, 3
          if(kdj(j,d) .eq. NANVAL) then
             ydj(j,d) = NANVAL
          else
             if(dd .eq. NANVAL) then
                ydj(j,d) = 0
             else
                ydj(j,d) = ydj(j,dd) + xj(kdj(j,dd))
             end if
             dd = d
          end if
       end do
       if(dd.ne.NANVAL) then
          xj(j) = ydj(j,dd) + xj(kdj(j,dd))
       end if
    end do

    uxj(1) = 1
    do j = 1, nj_
       dd = NANVAL
       do d = 0, 3
          if(ukdj(j,d) .eq. NANVAL) then
             uydj(j,d) = NANVAL
          else
             if(dd .eq. NANVAL) then
                uydj(j,d) = 0
             else
                uydj(j,d) = uydj(j, dd) + uxj(ukdj(j,dd))
             end if
             dd = d
          end if
       end do
       if(dd.ne.NANVAL) then
          uxj(j) = uydj(j,dd) + uxj(ukdj(j,dd))
       end if
    end do
    
  end subroutine DRTBuilder_xy
  function correct_abc(abc) result(res)
    integer, intent(in) :: abc(3)
    logical :: res
    res = all(abc.ge.0)
  end function correct_abc
end module Mod_DRTBuilder

!module Mod_DRTLoop
!  integer, num_seg_                           ! number of segments
!  integer, allocatable :: seg_table_(:,:,:)   ! (iseg, d,dd) gives segment-value
!  character(2), allocatable :: seg_names_(:)  ! (iseg) gives segment name
!  integer, parameter :: W=1, oR=2, uL=3, uR=4, oL=5, mR=6, pR=7, mL=8, pL=9
!  integer, parameter :: top_segs_(2) = (/oR,oL/)
!  integer, parameter :: bot_segs_(2) = (/uR,uL/)
!  integer, parameter :: mid_segs_(4) = (/mR,pR,mL,pL/)
!contains
!  subroutine DRTLoop_new
!    num_seg_ = 9
!    allocate(seg_table_(9,0:3,0:3))
!    allocate(seg_names_(num_seg_))
!  end subroutine DRTLoop_new
!end module Mod_DRTLoop

module Mod_DRT
  use Mod_ErrHandle
  use Mod_DRTConst, only : NANVAL
  implicit none
  type Obj_DRT
     integer :: nj, N, TwoS, norbs
     integer, allocatable :: abcj(:,:), kdj(:,:), xj(:), ydj(:,:)
     integer, allocatable :: ukdj(:,:), uxj(:), uydj(:,:)
  end type Obj_DRT
contains
  subroutine DRT_new(this, N, TwoS, norbs)
    use Mod_DRTBuilder
    type(Obj_DRT) :: this
    integer, intent(in) :: N, TwoS, norbs

    this%N = N
    this%TwoS = TwoS
    this%norbs = norbs
  
    ! -- determine structure of DRT --
    call DRTBuilder_new(N, TwoS, norbs)
    this%nj = nj_
    allocate(this%kdj(  this%nj, 0:3))
    allocate(this%abcj( this%nj, 3  ))
    allocate(this%xj(   this%nj     ))
    allocate(this%ydj(  this%nj, 0:3))
    allocate(this%ukdj( this%nj, 0:3))
    allocate(this%uxj(  this%nj     ))
    allocate(this%uydj( this%nj, 0:3))

    ! -- compute couting index and weight --
    call DRTBuilder_abck(this%abcj, this%kdj, this%ukdj)
    call DRTBuilder_xy(this%kdj, this%ukdj, &
         this%xj, this%ydj, this%uxj, this%uydj)

    ! -- finalize --
    call DRTBuilder_delete

  end subroutine DRT_new
  function num_wk(this) result(res)
    type(Obj_DRT) :: this
    integer res
    res = this%xj(1)
  end function num_wk
  function index_wk(this, ds) result(res)
    type(Obj_DRT) :: this
    integer, intent(in) :: ds(:)
    integer :: res
    integer :: i, j
    res = 1
    j = 1
    do i = this%norbs, 1, -1
       res = res + this%ydj(j, ds(i))
       j = this%kdj(j, ds(i))
    end do
  end function index_wk
!  subroutine DRT_last_wk(this, j0, ds)
!    ! Gives last walk whose tail is i in downwards order
!    type(Obj_DRT) :: this
!    integer, intent(in) :: j0
!    integer, intent(out) :: ds(:)
!    integer :: j, i, d
!    ds(:) = NANVAL
!    j = j0
!    do i = i0, norbs-1
!       do d = 0, 3
!          if(this%ukdj(j, d) .ne. NANVAL) then
!             ds(i+1) = d
!          end if
!       end do
!       if(ds(i+1) .eq. NANVAL) then
!          throw_err("error", 1)
!       end if
!    end do
!  end subroutine DRT_last_wk
!  subroutine DRT_first_uwk(this, j0, ds)
!    ! Gives first walk whose tail is i in upward order
!    type(Obj_DRT) :: this
!    integer, intent(in) :: j0
!    integer, intent(out) :: ds(:)
!    integer :: j, i, d
!    ds(:) = NANVAL
!    j = j0
!    do i = i0, this%norbs-1
!       do d = 3, 0, -1
!          if(this%kdj(j, d) .ne. NANVAL) then
!             ds(i) = d
!          end if
!       end do
!       if(ds(i) .eq. NANVAL) then
!          throw_err("error", 1)
!       end if
!    end do    
!  end subroutine DRT_first_uwk
!  subroutine get_wk(this, j_head, j_tail, d)
!    type(Obj_DRT) :: this
!    integer, intent(in) :: j_head, j_tail
!    integer, intent(out) :: d(this%norbs)
!    d(:) = NANVAL
!  end subroutine get_wk
  function ilevel(this, j) result(res)
    type(Obj_DRT) :: this
    integer, intent(in) :: j
    integer :: res
    res = sum(this%abcj(j,:))
  end function ilevel
  function correct_wk(this, j_head, j_tail, ds) result(res)
    type(Obj_DRT) :: this
    integer, intent(in) :: j_head, j_tail
    integer, intent(in) :: ds(this%norbs)
    integer :: i_head, i_tail, jj, i, j
    logical :: res
    
    res = .false.
    if(.not. j_head < j_tail) then
       throw_err("j_head < j_tail is required", 1)
    end if

    i_head = ilevel(this, j_head)
    i_tail = ilevel(this, j_tail)

    if(.not. i_head > i_tail) then
       throw_err("i_head < i_tail is required", 1)
    end if

    j = j_head
    do i = i_head, i_tail+1, -1
       if(ds(i) .eq. NANVAL) then
          throw_err("invalid Walk.", 1)
       end if
       jj = this%kdj(j, ds(i))
       if(jj.eq.NANVAL) then
          res = .false.
          return
       else
          j = jj
       end if
    end do

    if(j.eq.j_tail) then
       res = .true.
    else
       res = .false.
    end if
    
  end function correct_wk
  subroutine nodes_wk(this, ds, js)
    ! get node numbers {j} for walk {d}.
    type(Obj_DRT) :: this
    integer :: ds(this%norbs)
    integer :: js(0:this%norbs)
    integer j, i

    js(this%norbs) = 1
    do i = this%norbs, 1, -1
       j = this%kdj(js(i), ds(i))
       js(i-1) = j
    end do
    
  end subroutine nodes_wk
  subroutine nodes_sub_wk(this, j_head, ds, js)
    ! get node numbers {j} for walk {d}.
    type(Obj_DRT) :: this
    integer :: ds(this%norbs)
    integer :: j_head
    integer :: js(0:this%norbs)
    integer j, i, i_head

    i_head = ilevel(this, j_head)
    js(:) = NANVAL
    js(i_head) = j_head

    do i = i_head, 1, -1
       if(ds(i).eq.NANVAL) then
          js(i-1) = NANVAL
       else
          js(i-1) = this%kdj(js(i), ds(i))
       end if
    end do
    
  end subroutine nodes_sub_wk
  subroutine next_wk(this, ds, findq)
    type(Obj_DRT) :: this
    integer, intent(inout) :: ds(this%norbs)
    logical, intent(out) :: findq
    integer :: i, i0, d, jj
    integer :: js(0:this%norbs)

    findq = .false.

    ! -- determine current node number list --
    call get_nodes(this, ds, js)
    
    ! -- search most low level i which can increment d --
    i0 = NANVAL
    do i = 1, this%norbs
       do d = ds(i)+1, 3
          jj = this%kdj(js(i), d)
          if(jj.ne.NANVAL) then
             i0 = i
             ds(i0) = d
             js(i-1) = jj
             goto 111
          end if
       end do
    end do
    
111 continue
    if(i0.eq.NANVAL) then
       findq = .false.
       return
    end if

    ! -- from i0, search most low case numbers until end
    do i = i0-1, 1, -1
       do d = 0, 3
          jj = this%kdj(js(i), d)
          if(jj .ne. NANVAL) then
             ds(i) = d
             js(i-1) = jj
             exit
          end if
       end do
    end do

    findq = .true.
    return
    
  end subroutine next_wk
  subroutine next_sub_wk(this, j_head, i_tail, ds, findq)
    type(Obj_DRT) :: this
    integer, intent(in) :: j_head, i_tail
    integer, intent(inout) :: ds(0:this%norbs)
    logical, intent(out) :: findq

    integer :: i_head, i, j
    logical :: correct_wk
    
    i_head = sum(this%abcj(j_head,:))

    if(.not. i_head > i_tail+1) then
       throw_err("invalid j", 1)
    end if
    findq = .false.

    ! -- determine current node number list --
    call nodes_sub_wk(this, j_head, ds, js)
    
    ! -- search most low level i which can increment d --
    i0 = NANVAL
    do i = i_tail, i_head
       do d = ds(i)+1, 3
          jj = this%kdj(js(i), d)
          if(jj.ne.NANVAL) then
             i0 = i
             ds(i0) = d
             js(i-1) = jj
             goto 111
          end if
       end do
    end do
    
111 continue
    if(i0.eq.NANVAL) then
       findq = .false.
       return
    end if

    ! -- from i0, search most low case numbers until end
    do i = i0-1, i_tail, -1
       do d = 0, 3
          jj = this%kdj(js(i), d)
          if(jj .ne. NANVAL) then
             ds(i) = d
             js(i-1) = jj
             exit
          end if
       end do
    end do

    findq = .true.
    return
    
  end subroutine next_sub_wk
  !  subroutine DRT_loop1(this)
!    type(Obj_DRT) :: this
!    integer j_head, j_tail
!    integer :: ds_head(0:this%norbs), ds_tail(0:this%norbs)
!    integer :: ds_bra(0:this%norbs),  ds_ket(0:this%norbs)
!    integer :: dbra, dket, jbra, jket
!
!    do j_head = 1, this%nj
!       call DRT_last_wk(this, j_head, ds_head)
!       do j_tail = j_head, this%nj
!          if(sum(this%abcj(j_head,:)) - sum(this%abcj(j_tail,:)) < 2) then
!             cycle
!          end if
!          call DRT_first_wk(this, j_tail, ds_tail)
!          
!          ds_bra(:) = NANVAL
!          ds_ket(:) = NANVAL
!
!          jbra = j_head
!          jket = j_tail
!          do dbra = 0, 3
!             do dket = dbra+1, 3
!                if(this%kdj(jbra, dbra).ne.NANVAL .and. &
!                     this%kdj(jket, dket).ne.NANVAL ) then
!                end if
!             end do
!          end do
!          
!       end do
!    end do
!    
!  end subroutine DRT_loop1
  subroutine DRT_dump(this, downq, in_ifile)
    type(Obj_DRT) :: this
    logical, intent(in) :: downq
    integer, intent(in), optional :: in_ifile
    integer ifile, i, j
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    i = this%norbs
    write(ifile,*)  "  i |  j  |  aj  bj  cj | k0j k1j k2j k3j |  xj y0j y1j y2j y3j"
    write(ifile,*)  "----+-----+-------------+---------------------------------------"
    do j = 1, this%nj
       if(i .eq. sum(this%abcj(j,:))) then
          write(ifile,'("     |")', advance='no')
       else
          i = sum(this%abcj(j,:))
          write(ifile,'(i4," |")', advance='no') i
       end if

       write(ifile, '(i4, " |")', advance='no') j
       write(ifile, '(i4,i4,i4," |")', advance='no') this%abcj(j,:)
       if(downq) then
          write(ifile, '(i4,i4,i4,i4, " |")', advance='no') this%kdj(j,:)
          write(ifile, '(5i4)', advance='no') this%xj(j), this%ydj(j,:)
       else
          write(ifile, '(i4,i4,i4,i4, " |")', advance='no') this%ukdj(j,:)
          write(ifile, '(5i4)', advance='no') this%uxj(j), this%uydj(j,:)
       end if
       write(ifile, *)
    end do
    
  end subroutine DRT_dump
end module Mod_DRT
