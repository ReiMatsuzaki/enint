#include "macros.fpp"

#if defined REALGTO
#define FIELD double precision
#endif
#if defined HERMITEGTO
#define FIELD complex(kind(0d0)) 
#endif
#if defined COMPLEXGTO
#define FIELD complex(kind(0d0)) 
#endif
  
module Mod_Shel
  use Mod_ErrHandle
  implicit none
  type Obj_Shel
     FIELD, allocatable :: zeta(:)
     FIELD, allocatable :: coef(:,:) ! coef(mu,ig) 
     integer, allocatable :: ns(:,:) ! ns(i,:) = (nx,ny,nz) for i th basis
     integer :: ng, num, maxn
     FIELD :: w(3)
     integer :: ia ! index of atom
  end type Obj_Shel
contains
  subroutine Shel_new(this, ng, num)
    type(Obj_Shel) :: this
    integer, intent(in) :: ng, num

    this%ng = ng
    this%num = num
    allocate(this%zeta(ng))
    allocate(this%coef(num,ng))
    allocate(this%ns(num,3))

    this % zeta(:) = -1
    this % coef(:,:) = -1
    this % ns(:,:) = -1
    this % maxn = -1
    this % w(:) = 0
    
  end subroutine Shel_new
  subroutine Shel_ao_at(this, rs, n_dr, res)
    ! compute grid represented AO.
    ! j  : grid index
    ! mu : AO index
    type(Obj_Shel) :: this
    FIELD, intent(in) :: rs(:,:)
    integer, intent(in) :: n_dr(3)
    FIELD, intent(out) :: res(:,:)
    integer :: nj, n_mu, j, mu, ig, ir, ird
    FIELD, allocatable :: rj(:,:), r2j(:), exp_ig_j(:,:)
    FIELD :: tmp1, tmp2

    nj = size(rs, 1)
    n_mu = this%num
    allocate(rj(      nj, 3))
    allocate(r2j(     nj    ))

    do j = 1, nj
       rj(j,:) = rs(j,:) - this%w(:)
       r2j(j) = sum(rj(j,:)**2)
    end do

    allocate(exp_ig_j(this%ng, nj))
    do j = 1, nj
       do ig = 1, this%ng
          exp_ig_j(ig,j) = exp(-this%zeta(ig)*r2j(j))
       end do
    end do
    
    if(all(n_dr==0)) then
       do j = 1, nj
          do mu = 1, n_mu
             tmp1 = 1
             do ir = 1, 3
                tmp1 = tmp1 * rj(j,ir)**this%ns(mu,ir)
             end do
             tmp2 = 0
             do ig = 1, this%ng
                tmp2 = tmp2 + this%coef(mu,ig) * exp_ig_j(ig,j)
             end do
             res(mu,j) = tmp1*tmp2
          end do
       end do
    else if(sum(n_dr)==1) then
       ird = sum(n_dr*(/1,2,3/))
       do j = 1, nj
          do mu = 1, n_mu
             tmp1 = 1
             do ir = 1, 3
                tmp1 = tmp1 * rj(j,ir)**(this%ns(mu,ir) + n_dr(ir))
             end do
             tmp2 = 0
             do ig = 1, this%ng
                tmp2 = tmp2 -2*this%zeta(ig)* this%coef(mu,ig) * exp_ig_j(ig,j)
             end do
             res(mu,j) = tmp1*tmp2
             
             if(this%ns(mu,ird)>0) then
                tmp1 = this%ns(mu,ird)
                do ir = 1, 3
                   tmp1 = tmp1 * rj(j,ir)**(this%ns(mu,ir) - n_dr(ir))
                end do
                tmp2 = 0
                do ig = 1, this%ng
                   tmp2 = tmp2 + this%coef(mu,ig) * exp_ig_j(ig,j)
                end do
                res(mu,j) = res(mu,j) + tmp1*tmp2
             end if
          end do
       end do
    else
       throw_err("not impl", 1)
    end if
    
  end subroutine Shel_ao_at
  subroutine Shel_delete(this)
    type(Obj_Shel) :: this
    deallocate(this%zeta, this%coef, this%ns)
  end subroutine Shel_delete
end module Mod_Shel

module Mod_Nucs
  use Mod_ErrHandle
  implicit none
  type Obj_Nucs
     integer :: num
     FIELD, allocatable :: ws(:,:)
     FIELD, allocatable :: zs(:)
  end type Obj_Nucs
contains
  subroutine Nucs_new(this, num)
    type(Obj_Nucs) :: this
    integer, intent(in) :: num
    this%num = num
    allocate(this%ws(num,3))
    allocate(this%zs(num))
  end subroutine Nucs_new
  subroutine Nucs_delete(this)
    type(Obj_Nucs) :: this
    deallocate(this%ws, this%zs)
  end subroutine Nucs_delete
end module Mod_Nucs

Module Mod_Nshel
  use Mod_ErrHandle
  use Mod_Shel
  use Mod_Nucs
  use Mod_Molfunc
  implicit none
  type Obj_Nshel
     integer :: num
     type(Obj_Shel), allocatable :: shels(:)
     integer, allocatable :: j0s(:)
     type(Obj_Nucs) :: nucs
     integer :: nbasis
     logical :: setupq
  end type Obj_Nshel
contains
  subroutine Nshel_new(this, natom, num)
    type(Obj_Nshel) :: this
    integer, intent(in) :: num, natom
    this%num = num
    allocate(this%shels(num))
    allocate(this%j0s(num))
    call Nucs_new(this%nucs, natom); check_err()

    this%j0s(:) = -1
    this%setupq = .false.
    
  end subroutine Nshel_new
  subroutine Nshel_new_json(this, o)
    use Mod_fjson
    use Mod_Math
    type(Obj_Nshel) :: this
    type(object) :: o    
    FIELD, allocatable :: w(:,:), zs(:), coef_l(:,:)
    integer ng, num, natom, idx, ia, js, k0, k1
    integer, allocatable :: katom(:), kmin(:), kmax(:), kstart(:), kng(:)    
    character(5) :: ntypes(10)

    call object_get_i(o, "ng", ng); check_err()
    call object_get_i(o, "nshell", num); check_err()

    call object_get_idx(o, "katom", idx); check_err()
    allocate(katom(num))
    call a2ivec(o%vals(idx)%val_a, katom); check_err()
    natom = maxval(katom(:))

    ! -- allocate --
    call Nshel_new(this, natom, num); check_err()

    ! -- basic --
    allocate(kmin(num), kmax(num), kstart(num), kng(num))
    call object_get_idx(o, "kmin", idx); check_err()
    call a2ivec(o%vals(idx)%val_a, kmin(:)); check_err()
    call object_get_idx(o, "kmax", idx); check_err()
    call a2ivec(o%vals(idx)%val_a, kmax(:)); check_err()
    call object_get_idx(o, "kstart", idx); check_err()
    call a2ivec(o%vals(idx)%val_a, kstart(:)); check_err()
    call object_get_idx(o, "kng", idx); check_err()
    call a2ivec(o%vals(idx)%val_a, kng(:)); check_err()            
    
    ! -- nucs --
    call object_get_idx(o, "c", idx); check_err()
    allocate(w(3,natom))
    call a2dmat(o%vals(idx)%val_a, w(:,:)); check_err()
    this%nucs%ws(:,:) = transpose(w(:,:))

    call object_get_idx(o, "zan", idx); check_err()
    do ia = 1, natom
       this%nucs%zs(ia) = o%vals(idx)%val_a%vals(ia)%val_d
    end do

    
    ! -- GTO --
    allocate(zs(ng), coef_l(0:2,ng))
    
    call object_get_idx(o, "ex", idx); check_err()
    call a2dvec(o%vals(idx)%val_a, zs(:)); check_err()

    call object_get_idx(o, "cs", idx); check_err()
    call a2dvec(o%vals(idx)%val_a, coef_l(0,:)); check_err()
    call object_get_idx(o, "cp", idx); check_err()
    call a2dvec(o%vals(idx)%val_a, coef_l(1,:)); check_err()
    call object_get_idx(o, "cd", idx); check_err()
    call a2dvec(o%vals(idx)%val_a, coef_l(2,:)); check_err()
    
    do js = 1, num
       k0 = kstart(js)
       k1 = k0+kng(js)-1
       ntypes(:) = ""
       if(kmin(js)==1 .and. kmax(js)==1) then
          ntypes(1) = "s"
       else if(kmin(js)==1 .and. kmax(js)==4) then
          ntypes(1) = "s"; ntypes(2) = "p";
       else if(kmin(js)==2 .and. kmax(js)==4) then
          ntypes(1) = "p"
       else if(kmin(js)==5 .and. kmax(js)==10) then
          ntypes(1) = "d"
       else
          throw_err("unsupported combination", 1)
       end if
       call Nshel_set(this, js, ntypes, &
            kng(js), zs(k0:k1), coef_l(:,k0:k1), katom(js)) ; check_err()       
    end do
    
  end subroutine Nshel_new_json
  subroutine Nshel_new_file(this, fn)
    use Mod_fjson
    type(Obj_Nshel) :: this
    character(*), intent(in) :: fn
    type(value) :: v
    integer, parameter :: ifile = 12319291

    call loads_json_file(fn, ifile, v); check_err()
    close(ifile)
    call Nshel_new_json(this, v%val_o); check_err()
    call value_delete(v)
    
  end subroutine Nshel_new_file
  subroutine Nshel_delete(this)
    type(Obj_Nshel) :: this
    integer :: i
    do i = 1, this%num
       call Shel_delete(this%shels(i))
    end do
    deallocate(this%shels)
    call Nucs_delete(this%nucs)
  end subroutine Nshel_delete
  subroutine Nshel_dump(this, in_ifile)
    type(Obj_Nshel) :: this
    integer, intent(in) , optional :: in_ifile
    integer ifile, js, ia, j, jj, jg
    
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    write(ifile,*)
    write(ifile,*) "==== Nshel ===="
    write(ifile,*) "---- nucs  ----"
    write(ifile,'("ia",4A10)') "x", "y", "z", "Q"
    do ia = 1, this%nucs%num
       write(ifile,'(i0,4x,4f10.5)') ia, this%nucs%ws(ia,:), this%nucs%zs(ia)
    end do
    do js = 1, this%num
       write(ifile,*) 
       write(ifile,'(A,I0,A)') " ---- jshel = ", js, " ----"
       write(ifile, '("w    = ", 3f10.5)') this%shels(js)%w(:)
       write(ifile, '("zeta = ")', advance='no')
       do jg = 1, this%shels(js)%ng
          write(ifile,'(f10.5)',advance='no') this%shels(js)%zeta(jg)
       end do
       write(ifile, *)
       write(ifile,'(A4, A4, A10)') "j", "jj", "   n        coef"
       do jj = 1, this%shels(js)%num
          j = this%j0s(js)+jj
          write(ifile,'(I4,I4,3i3)',advance='no') j, jj, this%shels(js)%ns(jj,:)
          do jg = 1, this%shels(js)%ng
             write(ifile,'(2x,f10.5)',advance='no') this%shels(js)%coef(jj,jg)
          end do
          write(ifile,*) 
       end do
    end do
    
  end subroutine Nshel_dump
  subroutine Nshel_set(this, js, ntypes, ng, zeta, coef_l, ia)
    type(Obj_Nshel) :: this
    integer, intent(in) :: js, ng
    character(*), intent(in) :: ntypes(:)
    FIELD, intent(in) :: zeta(:)
    FIELD, intent(in) :: coef_l(0:,:)
    integer, intent(in) :: ia
    integer it, num, j, jj

    num = 0
    do it = 1, size(ntypes)
       if(ntypes(it) == "s") then
          num = num + 1
       else if(ntypes(it) == "p") then
          num = num + 3
       else if(ntypes(it) == "d") then
          num = num + 6
       else if(is_in_s(ntypes,"")) then
       else
          throw_err("unsupported", 1)
       end if
    end do

    call Shel_new(this%shels(js), ng, num); check_err()

    this%shels(js)%zeta(:) = zeta(1:ng)
    j = 1
    do it = 1, size(ntypes)
       if(ntypes(it) == "s") then
          this%shels(js)%ns(j,:) = (/0,0,0/)
          this%shels(js)%coef(j,:) = coef_l(0,1:ng)
          j = j + 1
       else if(ntypes(it) == "p") then
          this%shels(js)%ns(j+0,:) = (/1,0,0/)
          this%shels(js)%ns(j+1,:) = (/0,1,0/)
          this%shels(js)%ns(j+2,:) = (/0,0,1/)
          do jj = j, j+2
             this%shels(js)%coef(jj,:) = coef_l(1,1:ng)
          end do
          j = j + 3
       else if(ntypes(it) == "d") then
          this%shels(js)%ns(j+0,:) = (/2,0,0/)
          this%shels(js)%ns(j+1,:) = (/0,2,0/)
          this%shels(js)%ns(j+2,:) = (/0,0,2/)
          this%shels(js)%ns(j+3,:) = (/1,1,0/)
          this%shels(js)%ns(j+4,:) = (/1,0,1/)
          this%shels(js)%ns(j+5,:) = (/0,1,1/)
          do jj = j, j+5
             this%shels(js)%coef(jj,:) = coef_l(2,1:ng)
          end do
          j = j + 6
       else if(ntypes(it) == "") then
       else
          throw_err("unsupported", 1)
       end if
    end do
    
    this%shels(js)%maxn = maxval(this%shels(js)%ns(:,:))
    this%shels(js)%w(:) = this%nucs%ws(ia,:)
    this%shels(js)%ia = ia
    
  end subroutine Nshel_set
  subroutine Nshel_setup(this, in_normalize)
    type(Obj_Nshel) :: this
    logical ,intent(in), optional :: in_normalize
    integer :: j0, i, j, js, jj
    FIELD, allocatable :: smat(:,:)
    logical normalize
    double precision, parameter :: tol = 1.0d-10

    ! -- index --
    j0 = 0
    do i = 1, this%num
       this%j0s(i) = j0
       j0 = j0 + this%shels(i)%num
    end do
    this%nbasis = j0

    ! -- check current coefficient --
    do js = 1, this%num
       do jj = 1, this%shels(js)%num
          if(sum(abs(this%shels(js)%coef(jj,:)))<tol) then
             begin_err(1)
             write(0,*) "coefficient is too small"
             write(0,*) "jshel:", js
             write(*,*) "jj:", jj
             write(*,*) "sum(coef):", sum(abs( this%shels(js)%coef(jj,:)))
             end_err()
          end if
       end do
    end do

    ! -- normalization or not --
    if(present(in_normalize)) then
       normalize = in_normalize
    else
       normalize = .true.
    end if

    this%setupq = .true.
    
    if(normalize) then
       allocate(smat(this%nbasis, this%nbasis))
       call Nshel_s(this, smat); check_err()
       do js = 1, this%num
          do jj = 1, this%shels(js)%num
             j = this%j0s(js) + jj
             this%shels(js)%coef(jj,:) = this%shels(js)%coef(jj,:)/sqrt(smat(j,j))
          end do
       end do
    end if
    
  end subroutine Nshel_setup
  subroutine Nshel_s(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), i, j, k
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    FIELD :: d(3,0:5,0:5,0:10), acc, coef

    call check_setup(this); check_err()
    
    mat = 0
    do js = 1, this%num
       do ks = 1, this%num    
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk,0, d); check_err()

                do jj = 1, this%shels(js)%num
                   do kk = 1, this%shels(ks)%num                      
                      nj(:) = this%shels(js)%ns(jj,:)
                      nk(:) = this%shels(ks)%ns(kk,:)
                      acc = 1
                      do i = 1, 3
                         acc = acc * d(i,nj(i),nk(i),0)
                      end do
                      coef = cp * this%shels(js)%coef(jj,jg) &
                           * this%shels(ks)%coef(kk,kg)
                      j = this%j0s(js) + jj
                      k = this%j0s(ks) + kk
                      mat(j,k) = mat(j,k) + coef*acc
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine Nshel_s
  subroutine Nshel_t(this, mat)
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:)
    FIELD, allocatable :: mat_dw2(:,:,:)
    integer num

    call Nshel_num_basis(this, num)
    allocate(mat_dw2(3, num, num))
    call Nshel_dw2(this, mat_dw2)

    mat(:,:) = -0.5d0 * (mat_dw2(1,:,:) + mat_dw2(2,:,:) + mat_dw2(3,:,:))
    
  end subroutine Nshel_t
  subroutine Nshel_t_old(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), ir, jr, j, k, nkp(3)
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    FIELD :: acc, coef, cum    
    FIELD, allocatable :: d(:,:,:,:)

    call check_setup(this); check_err()

    allocate(d(1:3, 0:10, 0:10, 0:0))
    
    mat(:,:) = 0
    do js = 1, this%num
       do ks = 1, this%num
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk+2,0, d); check_err()

                do jj = 1, this%shels(js)%num
                   do kk = 1, this%shels(ks)%num
                      nj(:) = this%shels(js)%ns(jj,:)
                      nk(:) = this%shels(ks)%ns(kk,:)
                      
                      cum = 0
                      acc = 1
                      do ir = 1, 3
                         acc = acc * d(ir,nj(ir),nk(ir),0)
                      end do
                      cum = cum - 2*zk*(2*sum(nk)+3)*acc

                      do jr = 1, 3
                         nkp(:) = nk(:); nkp(jr) = nkp(jr)+2
                         acc = 1
                         do ir = 1, 3
                            acc = acc * d(ir,nj(ir),nkp(ir),0)
                         end do
                         cum = cum + 4*zk*zk*acc
                         
                         if(nk(jr)>1) then
                            nkp(:) = nk(:); nkp(jr) = nkp(jr)-2
                            acc = 1
                            do ir = 1, 3
                               acc = acc * d(ir,nj(ir),nkp(ir),0)
                            end do
                            cum = cum + nk(jr)*(nk(jr)-1)*acc
                         end if
                      end do
                      
                      coef = cp * this%shels(js)%coef(jj,jg) &
                           * this%shels(ks)%coef(kk,kg)
                      j = this%j0s(js) + jj
                      k = this%j0s(ks) + kk
                      mat(j,k) = mat(j,k) + coef*cum
                   end do
                end do
             end do
          end do
       end do
    end do
    mat(:,:) = -mat(:,:)/2
    deallocate(d)
  end subroutine Nshel_t_old
  subroutine Nshel_v(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), j, k, ic, maxn
    integer nx, ny, nz
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), wpc(3), ep, ccp, q
    FIELD :: d(3,0:5,0:5,0:10), acc, coef, cr(0:10,0:10,0:10)

    call check_setup(this); check_err()
    
    mat = 0
    do js = 1, this%num
       do ks = 1, this%num    
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                ccp = -2*pi*ep/zp
                call coef_d(zp,wp,wj,wk,maxnj,maxnk,maxnj+maxnk,d);
                check_err()

                do ic = 1, this%nucs%num
                   maxn = this%shels(js)%maxn + this%shels(ks)%maxn
                   wpc(:) = wp(:)-this%nucs%ws(ic,:)
                   q = this%nucs%zs(ic)
                   call coef_R(zp,wpc,maxn, cr); check_err()
                   do jj = 1, this%shels(js)%num
                      do kk = 1, this%shels(ks)%num                      
                         nj(:) = this%shels(js)%ns(jj,:)
                         nk(:) = this%shels(ks)%ns(kk,:)
                         
                         acc = 0
                         do nx = 0, nj(1)+nk(1)
                            do ny = 0, nj(2)+nk(2)
                               do nz = 0, nj(3)+nk(3)
                                  acc = acc + cr(nx,ny,nz) * q * &
                                       d(1,nj(1),nk(1),nx) * &
                                       d(2,nj(2),nk(2),ny) * &
                                       d(3,nj(3),nk(3),nz)
                               end do
                            end do
                         end do
                         coef = ccp * this%shels(js)%coef(jj,jg) &
                              * this%shels(ks)%coef(kk,kg)
                         j = this%j0s(js) + jj
                         k = this%j0s(ks) + kk
                         mat(j,k) = mat(j,k) + coef*acc
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine Nshel_v
  subroutine Nshel_r(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), j, k, ir, jr
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    FIELD :: d(3,0:5,0:5,0:10), coef, acc, acc0, acc1
    integer :: one(3,3)

    call check_setup(this); check_err()
    
    one = 0
    do ir = 1, 3
       one(ir,ir) = 1
    end do

    mat = 0
    do js = 1, this%num
       do ks = 1, this%num    
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk+1,0, d); check_err()

                do ir = 1, 3
                do jj = 1, this%shels(js)%num
                do kk = 1, this%shels(ks)%num                      
                   nj(:) = this%shels(js)%ns(jj,:)
                   nk(:) = this%shels(ks)%ns(kk,:)

                   acc0 = 1; acc1 = 1
                   do jr = 1, 3
                      acc0 = acc0 * d(ir,nj(jr),nk(jr),0)
                      acc1 = acc1 * d(ir,nj(jr),nk(jr)+one(ir,jr),0)
                   end do
                   acc = acc1 + wk(ir)*acc0
                   
                   coef = cp * this%shels(js)%coef(jj,jg) &
                        * this%shels(ks)%coef(kk,kg)
                   j = this%j0s(js) + jj
                   k = this%j0s(ks) + kk
                   mat(ir,j,k) = mat(ir,j,k) + coef*acc
                end do                   
                end do
                end do
             end do
          end do
       end do
    end do    
  end subroutine Nshel_r
  subroutine Nshel_dw(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), j, k, ir, jr
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    FIELD :: d(3,0:5,0:5,0:10), acc0, acc1, coef
    integer :: one(3,3)
    
    one = 0
    do ir = 1, 3
       one(ir,ir) = 1
    end do

    mat = 0
    do js = 1, this%num
       do ks = 1, this%num    
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk+1,0, d); check_err()

                do ir = 1, 3
                do jj = 1, this%shels(js)%num
                do kk = 1, this%shels(ks)%num                      
                   nj(:) = this%shels(js)%ns(jj,:)
                   nk(:) = this%shels(ks)%ns(kk,:)
                   acc0 = 2*zk
                   do jr = 1, 3
                      acc0 = acc0 * d(jr,nj(jr),nk(jr)+one(ir,jr),0)
                   end do
                   if(nk(ir)>0) then
                      acc1 = -nk(ir)
                      do jr = 1, 3
                         acc1 = acc1 * d(jr,nj(jr),nk(jr)-one(ir,jr),0)
                      end do
                   else
                      acc1 = 0
                   end if
                   
                   coef = cp* this%shels(js)%coef(jj,jg) * this%shels(ks)%coef(kk,kg)
                   j = this%j0s(js) + jj
                   k = this%j0s(ks) + kk
                   
                   mat(ir,j,k) = mat(ir,j,k) + coef*(acc0 + acc1)
                end do                   
                end do
                end do
             end do
          end do
       end do
    end do        
  end subroutine Nshel_dw
  subroutine Nshel_dw2(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: mat(:,:,:)
    integer num
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), ir, jr, j, k, nkp(3)  
    FIELD :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    FIELD :: acc, coef, cum    
    FIELD, allocatable :: d(:,:,:,:)

    call check_setup(this); check_err()

    call Nshel_num_basis(this, num)
    if(size(mat,1)<3 .or. size(mat,2)<num .or. size(mat,3)<num) then
       begin_err(1)
       write(0,*) "size mismatch"
       write(0,*) "basis size:", num
       write(0,*) "mat.shape:", shape(mat)
       end_err()
    end if

    allocate(d(1:3, 0:10, 0:10, 0:0))
    
    mat(:,:,:) = 0
    do js = 1, this%num
       do ks = 1, this%num
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk+2,0, d); check_err()

                do jj = 1, this%shels(js)%num
                   do kk = 1, this%shels(ks)%num
                      coef = cp * this%shels(js)%coef(jj,jg) &
                           * this%shels(ks)%coef(kk,kg)
                      j = this%j0s(js) + jj
                      k = this%j0s(ks) + kk
                      nj(:) = this%shels(js)%ns(jj,:)
                      nk(:) = this%shels(ks)%ns(kk,:)

                      do jr = 1, 3
                         cum = 0
                         
                         acc = 1
                         do ir = 1, 3
                            acc = acc * d(ir,nj(ir),nk(ir),0)
                         end do
                         cum = cum - 2*zk*(2*nk(jr)+1)*acc

                         nkp(:) = nk(:); nkp(jr) = nkp(jr)+2
                         acc = 1
                         do ir = 1, 3
                            acc = acc * d(ir,nj(ir),nkp(ir),0)
                         end do
                         cum = cum + 4*zk*zk*acc

                         if(nk(jr)>1) then
                            nkp(:) = nk(:); nkp(jr) = nkp(jr)-2
                            acc = 1
                            do ir = 1, 3
                               acc = acc * d(ir,nj(ir),nkp(ir),0)
                            end do
                            cum = cum + nk(jr)*(nk(jr)-1)*acc
                         end if
                         
                         mat(jr,j,k) = mat(jr,j,k) + coef*cum
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    deallocate(d)
  end subroutine Nshel_dw2
  subroutine Nshel_eri(this, eri)
    use Mod_Molfunc, only : coef_d
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    FIELD :: eri(:,:,:,:)
    integer is, js, ks, ls, ig, jg, kg, lg, ii, jj, kk, ll, i,j,k,l
    integer nxij, nyij, nzij, nxkl, nykl, nzkl, ni(3), nj(3), nk(3), nl(3)
    integer maxni, maxnj, maxnk, maxnl
    FIELD :: acc, cpp, coef, zi, zj, zk, zl, zij, zkl, eij, ekl
    FIELD :: wi(3), wj(3), wk(3), wl(3), wij(3), wkl(3), dij, dkl
    FIELD :: cdij(3,0:5,0:5,0:10), cdkl(3,0:5,0:5,0:10), cr(0:10,0:10,0:10)

    call check_setup(this); check_err()
    
    eri = 0
    do is = 1, this%num
    do js = 1, this%num
    do ks = 1, this%num
    do ls = 1, this%num
       maxni = this%shels(is)%maxn; wi = this%shels(is)%w(:)
       maxnj = this%shels(js)%maxn; wj = this%shels(js)%w(:)
       maxnk = this%shels(ks)%maxn; wk = this%shels(ks)%w(:)
       maxnl = this%shels(ls)%maxn; wl = this%shels(ls)%w(:)
       dij = dot_product(wi-wj,wi-wj);  dkl = dot_product(wk-wl,wk-wl); 
       do ig = 1, this%shels(is)%ng
       do jg = 1, this%shels(js)%ng
       do kg = 1, this%shels(ks)%ng
       do lg = 1, this%shels(ls)%ng
          zi = this%shels(is)%zeta(ig); zk = this%shels(ks)%zeta(kg)
          zj = this%shels(js)%zeta(jg); zl = this%shels(ls)%zeta(lg)
          zij = zi+zj; wij = (zi*wi+zj*wj)/zij; eij = exp(-zi*zj/zij*dij)
          zkl = zk+zl; wkl = (zk*wk+zl*wl)/zkl; ekl = exp(-zk*zl/zkl*dkl)
          cpp = 2*pi**(2.5d0)/(zij*zkl*sqrt(zij+zkl))*eij*ekl
          call coef_d(zij,wij,wi,wj,maxni,maxnj,maxni+maxnj, cdij)
          call coef_d(zkl,wkl,wk,wl,maxnk,maxnl,maxnk+maxnl, cdkl)
          call coef_R(zij*zkl/(zij+zkl),wi-wj,maxni+maxnj+maxnk+maxnl, cr)

          do ii = 1, this%shels(is)%num
          do jj = 1, this%shels(js)%num
          do kk = 1, this%shels(ks)%num
          do ll = 1, this%shels(ls)%num
             ni = this%shels(is)%ns(ii,:); nk = this%shels(ks)%ns(kk,:);
             nj = this%shels(js)%ns(jj,:); nl = this%shels(ls)%ns(ll,:);
             acc = 0
             do nxij = 0, ni(1)+nj(1)
             do nyij = 0, ni(2)+nj(2)
             do nzij = 0, ni(3)+nj(3)
             do nxkl = 0, nk(1)+nl(1)
             do nykl = 0, nk(2)+nl(2)
             do nzkl = 0, nk(3)+nl(3)
                acc = acc + cr(nxij+nxkl,nyij+nykl,nzij+nzkl) * &
                     cdij(1,ni(1),nj(1),nxij) * cdkl(1,nk(1),nl(1),nxkl) * &
                     cdij(2,ni(2),nj(2),nxij) * cdkl(2,nk(2),nl(2),nykl) * &
                     cdij(3,ni(3),nj(3),nxij) * cdkl(3,nk(3),nl(3),nzkl)
             end do
             end do
             end do
             end do
             end do
             end do
             coef = cpp * &
                  this%shels(is)%coef(ii,ig) * this%shels(ks)%coef(kk,kg) * &
                  this%shels(js)%coef(jj,jg) * this%shels(ls)%coef(ll,lg)
             i = this%j0s(is)+ii; k = this%j0s(ks)+kk
             j = this%j0s(js)+jj; l = this%j0s(ls)+ll
             eri(i,j,k,l) = eri(i,j,k,l) + coef*acc
             
          end do
          end do
          end do
          end do
                
       end do
       end do
       end do
       end do
    
    end do
    end do
    end do
    end do
 
  end subroutine Nshel_eri
  subroutine Nshel_vec_ia(this, vec_ia)
    type(Obj_Nshel) :: this
    integer, intent(out) :: vec_ia(:)
    integer :: num_basis, js ,j, jj

    call Nshel_num_basis(this, num_basis)
    if(size(vec_ia) < num_basis) then
       throw_err("not enough size", 1)
    end if

    do js = 1, this%num
       do jj = 1, this%shels(js)%num
          j = this%j0s(js) + jj
          vec_ia(j) = this%shels(js)%ia
       end do
    end do
    
  end subroutine Nshel_vec_ia
  subroutine Nshel_num_basis(this, num)
    type(Obj_Nshel) :: this
    integer, intent(out) :: num
    integer :: js
    
    num = 0
    do js = 1, this%num
       num = num + this%shels(js)%num
    end do
  end subroutine Nshel_num_basis
  subroutine check_setup(this)
    type(Obj_Nshel)::this
    if(.not. this%setupq) then
       throw_err("this object is not setup. call Nshel_setup first", 1)
    end if
  end subroutine check_setup
  subroutine Nshel_ao_at(this, rs, n_dr, res)
    type(Obj_Nshel) :: this
    double precision, intent(in) :: rs(:,:)
    integer, intent(in) :: n_dr(3)
    double precision, intent(out) :: res(:,:)

    integer mu0, mu1, ish
    
    call check_setup(this); check_err()

    mu0 = 1
    do ish = 1, this%num
       mu1 = mu0 + this%shels(ish)%num
       call Shel_ao_at(this%shels(ish), rs, n_dr, res(mu0:mu1-1,:)); check_err()
       mu0 = mu1
    end do
    
  end subroutine Nshel_ao_at
  ! == utils ==
  function is_in_s(ss, s) result(res)
    character(*), intent(in) :: ss(:)
    character(*) :: s
    logical :: res
    integer i
    res = .false.
    do i = 1, size(ss)
       if(trim(ss(i)) .eq. trim(s)) then
          res = .true.
       end if
    end do    
  end function is_in_s
end Module Mod_Nshel
