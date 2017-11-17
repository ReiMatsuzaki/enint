#include "macros.fpp"

module Mod_Nanval
  integer, parameter :: NANVAL=-1
end module Mod_Nanval

module Mod_DRTConst
  integer :: delta_abc_(0:3,3) = reshape((/0,0,1,1,  0, 1,-1, 0,  1, 0, 1, 0/), (/4,3/))
end module Mod_DRTConst

module Mod_DRTInitTable
  ! manage Distinct Row number j and (a,b,c) values.
  ! used for only initialization of DRT.
  use Mod_Nanval  
  implicit none
  integer, allocatable :: j_table_(:,:,:), abc_table_(:,:)
  integer :: current_j_
contains
  subroutine DRTInitTable_new(N, TwoS, norbs)
    integer, intent(in) :: N, TwoS, norbs
    integer :: na, nb, nc
    
    ! -- compute table --
    current_j_ = 1
    na = int(N/2)
    nb = int(TwoS + norbs)
    nc = int(norbs)

    allocate(j_table_(0:na, 0:nb, 0:nc))
    allocate(abc_table_((na+1)*(nb+1)*(nc+1),3))
    j_table_(:,:,:) = NANVAL
    abc_table_(:,:) = NANVAL
  end subroutine DRTInitTable_new
  subroutine DRTInitTable_delete
    deallocate(j_table_)
    deallocate(abc_table_)
  end subroutine DRTInitTable_delete
  subroutine get_j(abc, j)
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
    use Mod_DRTConst
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
!
!    if(correct_abc(aa,bb,cc)) then
!             call get_j(aa, bb, cc, jj)
!             this%kdj(j-j0+1,d) = jj
!          else
!             this%kdj(j-j0+1,d) = NANVAL
!          end if
    
  end subroutine next_j0j1
  function correct_abc(abc) result(res)
    integer, intent(in) :: abc(3)
    logical :: res
    res = all(abc.ge.0)
  end function correct_abc
end module Mod_DRTInitTable

module Mod_DRT
  use Mod_Nanval
  implicit none
  type Obj_DRT
     integer :: nj, N, TwoS, norbs
     integer, allocatable :: j0s(:), j1s(:)
     integer, allocatable :: abcj(:,:), kdj(:,:), xj(:), ydj(:,:)
  end type Obj_DRT
contains
  subroutine DRT_new(this, N, TwoS, norbs)
    use Mod_DRTInitTable
    type(Obj_DRT) :: this
    integer, intent(in) :: N, TwoS, norbs
    integer :: abc(3), i, j0, j1

    this%N = N
    this%TwoS = TwoS
    this%norbs = norbs
    allocate(this%j0s(0:norbs), this%j1s(0:norbs))

    ! -- determine structure of DRT using DRTInitTable --
    call DRTInitTable_new(N, TwoS, norbs)
    abc(1) = (N-TwoS)/2
    abc(2) = TwoS
    abc(3) = norbs-abc(1)-abc(2)
    call get_j(abc, this%j0s(norbs))
    this%j1s(norbs) = this%j0s(norbs)
    do i = norbs, 1, -1
       call next_j0j1(this%j0s(i), this%j1s(i), j0, j1)
       this%j0s(i-1) = j0; this%j1s(i-1) = j1
    end do

    ! -- allocation --
    this%nj = this%j1s(0)
    allocate(this%kdj(this%nj,0:3))
    allocate(this%abcj(this%nj,3))
    allocate(this%xj(this%nj))
    allocate(this%ydj(this%nj,0:3))

    ! -- compute (a,b,c) and (k_dj) --
    call new_abck(this)
    call DRTInitTable_delete

    ! -- compute couting indexes --
    call new_xy(this)

  end subroutine DRT_new
  subroutine new_abck(this)
    use Mod_DRTConst
    use Mod_DRTInitTable
    type(Obj_DRT) :: this
    integer :: j, d, abc_d(3), jj
    
    do j = 1, this%nj
       this%abcj(j,:) = abc_table_(j,:)
       do d = 0, 3
          abc_d(:) = this%abcj(j,:) - delta_abc_(d,:)
          if(correct_abc(abc_d)) then
             call get_j(abc_d, jj)
             this%kdj(j,d) = jj
          else
             this%kdj(j,d) = NANVAL
          end if
       end do
    end do
    
  end subroutine new_abck
  subroutine new_xy(this)
    type(Obj_DRT) :: this

    integer j, dd, d
    this%xj(this%nj) = 1
    do j = this%nj, 1, -1
       dd = NANVAL
       do d = 0, 3
          if(this%kdj(j,d) .eq. NANVAL) then
             this%ydj(j,d) = NANVAL
          else
             if(dd .eq. NANVAL) then
                this%ydj(j,d) = 0
             else
                this%ydj(j,d) = this%ydj(j,dd) + this%xj(this%kdj(j,dd))
             end if
             dd = d
          end if
       end do
       if(dd.ne.NANVAL) then
          this%xj(j) = this%ydj(j,dd) + this%xj(this%kdj(j,dd))
       end if
    end do
    
  end subroutine new_xy
  subroutine DRT_dump(this, in_ifile)
    type(Obj_DRT) :: this
    integer, intent(in), optional :: in_ifile
    integer ifile, i, k0, k1, k2, k3, j
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    write(ifile,*)  "  i |  j  |  aj  bj  cj | k0j k1j k2j k3j |  xj y0j y1j y2j y3j"
    write(ifile,*)  "----+-----+-------------+---------------------------------------"
    do i = this%norbs, 0, -1
       do j = this%j0s(i), this%j1s(i)
          k0 = this%kdj(j,0)
          k1 = this%kdj(j,1)
          k2 = this%kdj(j,2)
          k3 = this%kdj(j,3)
          if(j.eq.this%j0s(i)) then
             write(ifile,'(i4," |")', advance='no') i
          else
             write(ifile,'("     |")', advance='no') 
          end if
          write(ifile, '(i4, " |")', advance='no') j
          write(ifile, '(i4,i4,i4," |")', advance='no') this%abcj(j,:)
          write(ifile, '(i4,i4,i4,i4, " |")', advance='no') this%kdj(j,:)
          write(ifile, '(5i4)', advance='no') this%xj(j), this%ydj(j,:)
          write(ifile, *)
       end do
    end do
    
  end subroutine DRT_dump
end module Mod_DRT

!module Mod_DRT
!  use Mod_PreDRT
!  implicit none
!contains
!  subroutine DRT_new()
!end module Mod_DRT

module Mod_PreDRT_old
  use Mod_ErrHandle
  use Mod_Nanval
  implicit none
  type Obj_DRs
     integer :: num
     integer, allocatable :: jabcs(:)
     integer, allocatable :: kdj(:,:)
  end type Obj_DRs
  integer :: norbs_
  integer :: delta_a(0:3) = (/0, 0, 1, 1/)
  integer :: delta_b(0:3) = (/0, 1,-1, 0/)
  integer :: delta_c(0:3) = (/1, 0, 1, 0/)
  integer, allocatable :: jabc_table_(:,:,:)
  integer, allocatable :: abc_table_(:,:)
  type(Obj_DRs), allocatable :: drs_(:)
contains
  subroutine PreDRT_new(N, TwoS, norbs)
    integer, intent(in) :: N, TwoS
    integer, intent(in) :: norbs
    integer :: a, b, c, i

    norbs_ = norbs

    call make_table(N, TwoS, norbs)
    allocate(drs_(0:norbs))

    a = (N-TwoS)/2
    b = TwoS
    c = norbs-a-b

    call allocate_dr(norbs, 1); check_err()
    drs_(norbs)%jabcs(1) = jabc_table_(a,b,c)
    call add_k(norbs)

    do i = norbs-1, 0, -1
       call add_jabc(i)
       call add_k(i)
    end do
    
  end subroutine PreDRT_new
  subroutine PreDRT_delete
    integer i
    do i = 0, norbs_
       deallocate(drs_(i)%jabcs)
       deallocate(drs_(i)%kdj)
    end do
    deallocate(drs_, jabc_table_, abc_table_)
    
  end subroutine PreDRT_delete
  subroutine make_table(N, TwoS, norbs)
    integer, intent(in) :: N, TwoS, norbs
    integer :: na, nb, nc
    integer :: jabc, a, b, c
    
    ! -- compute table --
    na = int(N/2)
    nb = int(TwoS + norbs)
    nc = int(norbs)

    allocate(jabc_table_(0:na, 0:nb, 0:nc))
    allocate(abc_table_((na+1)*(nb+1)*(nc+1),3))

    jabc = 1
    do a = 0, na
       do b = 0, nb
          do c = 0, nc
             jabc_table_(a,b,c) = jabc
             abc_table_(jabc,1) = a
             abc_table_(jabc,2) = b
             abc_table_(jabc,3) = c
             jabc = jabc + 1
          end do
       end do
    end do

  end subroutine make_table
  subroutine allocate_dr(i, ndr)
    integer, intent(in) :: i
    integer, intent(in) :: ndr
    drs_(i)%num = ndr
    allocate(drs_(i)%jabcs(ndr))
    allocate(drs_(i)%kdj(ndr,0:3))
    drs_(i)%jabcs = NANVAL
    drs_(i)%kdj = NANVAL
  end subroutine allocate_dr
  function correct_abc(a,b,c) result(res)
    integer, intent(in) :: a, b, c
    logical :: res
    res = (a.ge.0 .and. b.ge.0 .and. c.ge.0)
  end function correct_abc
  subroutine add_k(i)
    integer, intent(in) :: i
    integer :: jabc, idr, d, aj, bj, cj, jabc_new
    
    do idr = 1, drs_(i)%num
       jabc = drs_(i)%jabcs(idr)
       do d = 0, 3
          aj = abc_table_(jabc, 1) - delta_a(d)
          bj = abc_table_(jabc, 2) - delta_b(d)
          cj = abc_table_(jabc, 3) - delta_c(d)
          if(correct_abc(aj,bj,cj)) then             
             jabc_new = jabc_table_(aj,bj,cj)
             drs_(i)%kdj(idr,d) = jabc_new
          else
             drs_(i)%kdj(idr,d) = NANVAL
          end if
       end do
    end do
    
  end subroutine add_k
  subroutine add_jabc(i)
    integer, intent(in) :: i
    integer :: idx, d, num_dr_i, idxi, jabc

    num_dr_i = 0
    do idx = 1, drs_(i+1)%num
       do d = 0, 3
          jabc = drs_(i+1)%kdj(idx,d)
          if(jabc .ne. NANVAL) then             
             num_dr_i = num_dr_i + 1
          end if
       end do
    end do
    
    call allocate_dr(i, num_dr_i)
    
    idxi = 1
    do idx = 1, drs_(i+1)%num
       do d = 0, 3
          jabc = drs_(i+1)%kdj(idx,d)
          if(jabc .ne. NANVAL) then             
             drs_(i)%jabcs(idxi) = drs_(i+1)%kdj(idx,d)
             idxi = idxi+1
          end if
       end do
    end do
    
  end subroutine add_jabc
  subroutine PreDRT_dump(in_ifile)
    integer, intent(in), optional :: in_ifile
    integer ifile, i, idx, jabc, a, b, c, k0, k1, k2, k3
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    write(ifile,*)  "  i | jabc|  aj  bj  cj | k0j k1j k2j k3j"
    write(ifile,*)  "----+-----+-------------+------------------"
    do i = norbs_, 0, -1
       do idx = 1, drs_(i)%num
          jabc = drs_(i)%jabcs(idx)
          a = abc_table_(jabc,1)
          b = abc_table_(jabc,2)
          c = abc_table_(jabc,3)
          k0 = drs_(i)%kdj(idx,0)
          k1 = drs_(i)%kdj(idx,1)
          k2 = drs_(i)%kdj(idx,2)
          k3 = drs_(i)%kdj(idx,3)
          if(idx.eq.1) then
             write(ifile,'(i4," |")', advance='no') i
          else
             write(ifile,'("     |")', advance='no') 
          end if
          write(ifile, '(i4, " |")', advance='no') jabc
          write(ifile, '(i4,i4,i4," |")', advance='no') a,b,c
          write(ifile, '(i4,i4,i4,i4)', advance='no') k0,k1,k2,k3
          write(ifile, *)
       end do
    end do
    
  end subroutine PreDRT_dump
end module Mod_PreDRT_Old

module Mod_DRT_old
  use Mod_Nanval
  use Mod_ErrHandle  
  implicit none
  type Obj_DRT
     integer :: N, TwoS, norbs
     integer :: nj
     integer, allocatable :: j0(:), j1(:)
     integer, allocatable :: abcj(:,:), xj(:), ydj(:,:), kdj(:,:)
  end type Obj_DRT
contains
  subroutine DRT_new(this, N, TwoS, norbs)
    
    type(Obj_DRT) :: this
    integer, intent(in) :: N, TwoS, norbs

    call new_graph(this, N, TwoS, norbs)
    call new_weight(this)
    
  end subroutine DRT_new
  subroutine new_graph(this, N, TwoS, norbs)
    use Mod_PreDRT_old
    type(Obj_DRT) :: this
    integer, intent(in) :: N, TwoS, norbs
    integer :: i,j,ndr,d,idx,jabc
    integer, allocatable :: jtable(:)

    this%norbs = norbs
    call PreDRT_new(N, TwoS, norbs)

    this%nj = 0
    do i = 0, norbs
       this%nj = this%nj + drs_(i)%num
    end do

    allocate(this % j0(0:norbs))
    allocate(this % j1(0:norbs))
    allocate(this % abcj(this%nj,3))
    allocate(this % xj(this%nj))
    allocate(this % ydj(this%nj,0:3))
    allocate(this % kdj(this%nj,0:3))

    allocate(jtable(size(abc_table_, 1)))
    j = 1
    do i = norbs, 0, -1
       do idx = 1, drs_(i)%num
          jabc = drs_(i)%jabcs(idx)
          jtable(jabc) = j
       end do
    end do

    j = 1
    this%j0(norbs) = j
    do i = norbs, 1, -1
       ndr = drs_(i)%num
       this%j1(i) = this%j0(i)+ndr-1
       this%j0(i-1) = this%j0(i)+ndr
       do idx = 1, ndr
          jabc = drs_(i)%jabcs(idx)
          this%abcj(j,:) = abc_table_(jabc,:)
          do d = 0, 3
             if(drs_(i)%kdj(idx,d) .ne. NANVAL) then
                this%kdj(j,d) = jtable(drs_(i)%kdj(idx,d))
             else
                this%kdj(j,d) = NANVAL
             end if
          end do
          j = j + 1
       end do
    end do
    this%j1(0) = this%j0(0)
    this%abcj(j,:) = 0
    this%kdj(j,:) = NANVAL
    
  end subroutine new_graph
  subroutine new_weight(this)
    type(Obj_DRT) :: this
    integer i, j, dd, d
    
    this%xj(this%nj) = 1
    do i = 1, this%norbs       
       do j = this%j0(i), this%j1(i)
          dd = NANVAL
          do d = 0, 3
             if(this%kdj(j,d) .eq. NANVAL) then
                this%ydj(j,d) = NANVAL
             else
                if(d .eq. NANVAL) then
                   this%ydj(j,d) = 0
                else
                   this%ydj(j,d) = this%ydj(j,dd) + this%xj(this%kdj(j,dd))
                end if
                dd = d
             end if
          end do
          if(dd.ne.NANVAL) then
             this%xj(j) = this%ydj(j,dd) + this%xj(this%kdj(j,dd))
          end if
       end do
    end do
    
  end subroutine new_weight
  subroutine DRT_dump(this, in_ifile)
    type(Obj_DRT) :: this
    integer, intent(in), optional :: in_ifile
    integer ifile, i, j, a, b, c, k0, k1, k2, k3
    
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    write(ifile,*)  "  i |   j |  aj  bj  cj | k0j k1j k2j k3j"
    write(ifile,*)  "----+-----+-------------+------------------"
    do i = this%norbs, 0, -1
       do j = this%j0(i), this%j1(i)
          a = this%abcj(j,1)
          b = this%abcj(j,2)
          c = this%abcj(j,3)
          k0 = this%kdj(j,0)
          k1 = this%kdj(j,1)
          k2 = this%kdj(j,2)
          k3 = this%kdj(j,3)
          if(j.eq.this%j0(i)) then
             write(ifile,'(i4," |")', advance='no') i
          else
             write(ifile,'("     |")', advance='no') 
          end if
          write(ifile, '(i4, " |")', advance='no') j
          write(ifile, '(i4,i4,i4," |")', advance='no') a,b,c
          write(ifile, '(i4,i4,i4,i4)', advance='no') k0,k1,k2,k3
          write(ifile, *)
       end do
    end do
        
  end subroutine DRT_dump
end module Mod_DRT_Old
