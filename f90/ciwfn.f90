#include "macros.fpp"

module Mod_CIWfn
  use Mod_ErrHandle
  implicit none
  integer :: nloop_, nijIJ_, nfrozen_, nmo_, ncsf_
  integer, allocatable :: kstart_(:), kend_(:)
  integer, allocatable :: i_(:), j_(:), ii_(:), jj_(:)
  double precision, allocatable :: v_(:)
  end type Obj_Aij
  subroutine Aij_new(nloop, nijIJ, nfrozen)
    integer, intent(in) :: nloop, nijIJ, nfrozen
    nloop_ = nloop
    nijIJ_ = nijIJ
    nfrozen_ = nfrozen
    allocate(kstart_(nloop))
    allocate(kend_(nloop))
    allocate(kend_(nloop))
    allocate(i_(nijIJ))
    allocate(j_(nijIJ))
    allocate(ii_(nijIJ))
    allocate(jj_(nijIJ))
    allocate(v_(nloop))
  end subroutine Aij_new
  subroutine Aij_new_csv(ifile)
    integer, intent(in) :: ifile
    integer :: ijIJ, kloop, i, j, ii, jj, minj
    double precision :: val0, val
    double precision, parameter :: tol = 1.0d-5

    val0 = 0.0d0
    ijIJ = 0
    kloop = 0
    minj = -1
    read(ifile,*)
    do
       read(ifile,*, end=100) i,j,ii,jj,val
       
       if(minj.eq.-1 .or. nfrozen>i) then
          minj = i
       end if

       ijIJ = ijIJ + 1
       if(abs(val-val0) > tol) then
          kloop = kloop + 1
          val0 = val
       end if
    end do

100 continue    
    call Aij_new(kloop, ijIJ, minj-1); check_err()
    val0 = 0.0d0
    kstart_(1) = 1
    kend_(nloop_) = nijIJ_
    ijIJ = 0
    kloop = 0
    rewind(ifile)
    read(ifile,*)
    do
       read(ifile,*, end=200) i,j,ii,jj,val
       ijIJ = ijIJ + 1
       i_(ijIJ) = i
       j_(ijIJ) = j
       ii_(ijIJ) = ii
       jj_(ijIJ) = jj
       if(abs(val-val0) > tol) then
          kloop = kloop + 1
          val0 = val
          if(kloop.ne.1) then
             kend_(kloop-1) = ijIJ
          end if
          if(kloop.ne.nloop) then
             kstart_(kloop) = ijIJ
          end if
          v_(kloop) = val
       end if
    end do

200 continue
    
  end subroutine Aij_new_csv
  subroutine Aij_setup
    nmo_ = maxval(i_(:))
    ncsf_ = maxval(ii_(:))
  end subroutine Aij_setup
  subroutine Aij_dm1(cci, dm1)
    complex(kind(0d0)), intent(in) :: cci(:)
    complex(kind(0d0)), intent(out) :: dm1(:,:)
    integer ijIJ, kloop, i, j, ii, jj, ijIJ
    double precision v

    if(size(dm1,1).ne.nmo_ .or. size(dm1,2).ne.nmo_) then
       throw_err("invalid size",1)
    end if

    dm1 = 0
    do i = 1, nfrozen_
       do ii = 1, ncsf_
          dm1(i,i) = dm1(i,i) + 2*conjg(cci(ii))*cci(ii)
       end do
    end do
    
    do kloop = 1, nloop_
       do ijIJ = kstart_(kloop), kend_(kloop)
          i = i_( ijIJ)
          j = j_( ijIJ)
          ii= ii_(ijIJ)
          jj= jj_(ijIJ)
          v = v_( ijIJ)
          if(ii.eq.jj) then
             dm1(i,i) = dm1(i,i) + 2*conjg(cci(ii))*cci(ii)
          else
             dm1(i,j) = dm1(i,j) + v*conjg(cci(ii))*cci(jj)
             dm1(j,i) = dm1(j,i) + v*conjg(cci(jj))*cci(ii)
          end if
       end do
    end do
    
  end subroutine Aij_dm1
  subroutine ciwfn_op1(mmo, dm1, val)
    complex(kind(0d0)), intent(in) :: mmo(:,:)
    complex(kind(0d0)), intent(in) :: dm1(:,:)
    complex(kind(0d0)), intent(out) :: val

    val = sum(mmo(1:nmo_,1:nmo_)*dm1(1:nmo_,1:nmo_))
    
  end subroutine ciwfn_op1
  subroutine ao2mo(mao, cmo, mmo)
    complex(kind(0d0)), intent(in) :: mao(:,:), cmo(:,:)
    complex(kind(0d0)), intent(out) :: mmo(:,:)

    mmo = dot_product(conjg(transpose(cmo)), dot_product(mao, cmo))
    
  end subroutine ao2mo
  subroutine rho_mo(psi_ao, cmo, res)
    complex(kind(0d0)), intent(in) :: psi_ao(:)
    complex(kind(0d0)), intent(in) :: cmo(:,:)
    complex(kind(0d0)), intent(out) :: res(:,:,:)
    complex(kind(0d0)) :: psi_mo(size(cmo(0,:)))
    integer :: i, j, r

    psi_mo = dot_product(transpose(cmo), psi_ao)

    res = 0
    do i = 1, nmo_
       do j = 1, nmo_
          res(i,j) += conjg(psi_mo(i,:)) * psi_mo(j,:)
       end do
    end do
    
  end subroutine rho_mo
contains
  
end module Mod_CIWfn
