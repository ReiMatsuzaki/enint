#include "macros.fpp"

#if defined REALFUNC
#define FIELD double precision
#endif
#if defined COMPLEXFUNC
#define FIELD complex(kind(0d0)) 
#endif

module Mod_Molfunc
  use Mod_ErrHandle
  implicit none
  integer :: method_coef_R_ = 1
  integer :: method_igamma_ = 0
contains
  ! ==== Incomplete Gamma Fucntion ====
  subroutine igamma_c_f1(maxm, z, res)
    use Mod_Const, only : ii, pi
    integer, intent(in) :: maxm
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(0:)
    double precision, parameter :: eps = 1.0d-10
    double precision, parameter :: delta = 1.0d-15
    integer, parameter :: NF_MAX=50, NF_0=10
    integer :: m, n
    double precision :: anR(0:NF_MAX), anI(0:NF_MAX), z2, az, c, s, x, y
    double precision :: bmR(0:maxm), bmI(0:maxm), fmR(0:maxm), fmI(0:maxm)
    logical :: convq

    x = real(z)
    y = aimag(z)
    z2 = x*x+y*y
    az = abs(z)
    c = sqrt(pi*(az+x)/(8*z2))
    s = sqrt(pi*(az-x)/(8*z2)) + eps

    if(x < eps .or. y < -eps) then
       throw_err("Re[z] and Im[z] must be positive", 1)
    end if

    anR(:) = 0
    anI(:) = 0
    anR(0) = -x/(2*z2)*exp(-x)*cos(y) + y/(2*z2)*exp(-x)*sin(y)
    anI(0) = +y/(2*z2)*exp(-x)*cos(y) + x/(2*z2)*exp(-x)*sin(y)
    do n = 1, NF_MAX
       anR(n) = -(2*n-1) * (x/(2*z2)*anR(n-1) + y/(2*z2)*anI(n-1))
       anI(n) = +(2*n-1) * (y/(2*z2)*anR(n-1) - x/(2*z2)*anI(n-1))
       if(n>NF_0) then
          if(c*delta > abs(anR(n)) .and. s*delta>abs(anI(n))) then
             convq = .True.
             exit
          end if
       end if
    end do
    if(.not. convq) then
       throw_err("not converged",1)
    end if
    
    fmR(0) = +c + sum(anR)
    fmI(0) = -s + sum(anI)
    do m = 1, maxm
       fmR(m) = (2*m-1) * (x/(2*z2)*fmR(m-1) + y/(2*z2)*fmI(m-1))
       fmI(m) = (2*m-1) * (x/(2*z2)*fmI(m-1) - y/(2*z2)*fmR(m-1))
    end do

    bmR(0) = 0
    bmI(0) = 0
    do m = 1, maxm
       bmR(m) = anR(0) + (2*m-1)*(x/(2*z2)*bmR(m-1) + y/(2*z2)*bmI(m-1))
       bmI(m) = anI(0) + (2*m-1)*(x/(2*z2)*bmI(m-1) - y/(2*z2)*bmR(m-1))
    end do

    res(:) = (fmR(:) + bmR(:)) + ii * (fmI(:)+bmI(:))
    
  end subroutine igamma_c_f1
  recursive subroutine igamma_c_f2(maxm, z, res)
    integer, intent(in) :: maxm
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(0:)
    double precision, parameter :: eps = 1.0d-10
    integer, parameter :: NR=47
    integer :: m, n
    complex(kind(0d0)) :: bn(0:NR), an(0:NR)
    double precision tmp1, tmp3, tmp5, tmp7, e, f1, f2, f3, t1, t2, t3, x, y

    x = real(z)
    y = aimag(z)

    if(x < -eps) then
       throw_err("Re[z] must be positive", 1)
    end if
    if(y < -eps) then
       call igamma_c_f2(maxm, z, res); check_err()
       res(:) = conjg(res(:))
       return 
    end if

    bn(0) = 1
    bn(1) = 1 + z/2
    do n = 2, NR
       bn(n) = bn(n-1) + z*z/(4*(2*n-1)*(2*n-3))*bn(n-2)
    end do

    do m = 0, maxm
       an(0) = 1
       tmp1 = 2*m+1; tmp3 = tmp1+2; tmp5 = tmp3+2; tmp7 = tmp5+2       
       t1 = tmp1/tmp3
       t2 = tmp1/(tmp3*tmp5)
       t3 = tmp1*(tmp1**2+44)/(60*tmp3*tmp5*tmp7)
       an(1) = bn(1) -t1*z
       an(2) = bn(2) -t1*z -t2*z**2
       an(3) = bn(3) -t1*z -t2*z**2  -t3*z**3
       do n = 4, NR
          f1 = (2.0d0*n-2*m-5)/(2*(2*n-3)*(2*n+2*m+1))
          f2 = 1.0d0/(4*(2*n-1)*(2*n-3))
          f3 = -f1/(4*(2*n-3)*(2*n-5))
          e = -f1
          an(n) = (1+f1*z)*an(n-1) + (e+f2*z)*z*an(n-2) + f3*z**3*an(n-3)
       end do
       res(m) = 1.0d0/(2*m+1) * an(NR)/bn(NR)
    end do
    
  end subroutine igamma_c_f2
  recursive subroutine igamma_c(maxm, z, res)
    integer, intent(in) :: maxm
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(0:)
    double precision :: x, y
    double precision, parameter :: eps = 1.0d-15
    
    x = real(z)
    y = aimag(z)
    if(x > -eps .and. y > -eps) then
       if(x < 21 .and. x+y < 37) then
          call igamma_c_f2(maxm, z, res); check_err()
       else
          call igamma_c_f1(maxm, z, res); check_err()
       end if
    else
       call igamma_c(maxm, conjg(z), res(:)); check_err()
       res(:) = conjg(res(:))
    end if
    
  end subroutine igamma_c
  recursive function igamma_r_1(m, z) result(res)
    ! compute incompute gamma function 
    !        F_m(z) = Int_0^1 t^{2m} Exp[-zt^2] dt
    ! implemented function P in gsl is defined as
    !        P(a,x) = 1/Gamma(a) . Int_0^x t^{a-1}t^{-t} dt
    use fgsl
    integer, intent(in) :: m
    FIELD, intent(in) :: z
    FIELD :: res
    FIELD :: a
    real(fgsl_double) :: ga, giaz
    double precision, parameter :: eps=10.0d-14

    if(abs(z)<eps) then
       res = 1/(2*m+1.0d0)
       return
    end if

    a = m+0.5d0
    ga = fgsl_sf_gamma(a)
    giaz = fgsl_sf_gamma_inc_P(a, z)
    res = ga/(2*z**a) * giaz
    return

  end function igamma_r_1
  subroutine igamma_r_fast(maxm, z, res)
    use Mod_const, only : pi
    integer, intent(in) :: maxm
    FIELD, intent(in) :: z
    FIELD, intent(out) :: res(0:)
    integer m

    if(abs(z)<1.0d-4) then
       do m = 0, maxm
          res(m) = 1/(2*m+1.0d0) &
               - z/(2*m+3) &
               + (z**2)/(2*(2*m+5)) &
               - (z**3)/(3*(2*m+7))
       end do
    else
       res(0) = sqrt(pi/z)/2 * erf(sqrt(z))
       do m = 1, maxm
          res(m) = -exp(-z)/(2*z) + (2*m-1)/(2*z) * res(m-1)
       end do
    end if
    
  end subroutine igamma_r_fast
  subroutine igamma_r(maxm, z, res)
    integer, intent(in) :: maxm
    FIELD, intent(in) :: z
    FIELD :: res(0:);
    integer m

    select case(method_igamma_)
    case(0)
       do m = 0, maxm
          res(m) = igamma_r_1(m, z)
       end do
    case(1)
       call igamma_r_fast(maxm, z, res); check_err()
    end select
    
  end subroutine igamma_r
  subroutine igamma(maxm, z, res)
    integer, intent(in) :: maxm
    FIELD, intent(in) :: z
    FIELD, intent(out) :: res(0:)

#if defined REALFUNC
    call igamma_r(maxm, z, res); check_err()
#endif    
#if defined COMPLEXFUNC
    call igamma_c(maxm, z, res); check_err()    
#endif    
  end subroutine igamma
  ! ==== Hermitian coefficient ====
  recursive function coef_d1(zp,wp,wj,wk,nj,nk,n) result(res)
    FIELD, intent(in) :: zp, wp, wj, wk
    integer, intent(in) :: nj, nk, n
    FIELD :: res

    if(nj==0 .and. nk==0 .and. n==0) then
       res = 1.0
    else if(n<0 .or. n>nj+nk) then
       res = 0.0
    else if(nj>0) then
       res = 1/(2*zp) * coef_d1(zp,wp,wj,wk,nj-1,nk,n-1) + &
            (wp-wj)   * coef_d1(zp,wp,wj,wk,nj-1,nk,n)   + &
            (n+1)     * coef_d1(zp,wp,wj,wk,nj-1,nk,n+1)
    else
       res = 1/(2*zp) * coef_d1(zp,wp,wj,wk,nj,nk-1,n-1) + &
            (wp-wk)   * coef_d1(zp,wp,wj,wk,nj,nk-1,n)   + &
            (n+1)     * coef_d1(zp,wp,wj,wk,nj,nk-1,n+1)
    end if
    
  end function coef_d1
  subroutine coef_d(zp,wp,wj,wk,maxnj,maxnk,maxn,res)
    FIELD, intent(in) :: zp, wp(3), wj(3), wk(3)
    integer, intent(in) :: maxnj, maxnk, maxn
    FIELD, intent(out) :: res(:,0:,0:,0:)
    integer i, nj, nk, n

    do i = 1, 3
       do nj = 0, maxnj
          do nk = 0, maxnk
             do n = 0, maxn
                res(i,nj,nk,n) = coef_d1(zp,wp(i),wj(i),wk(i),nj,nk,n)
             end do
          end do
       end do
    end do
    
  end subroutine coef_d
  ! ==== Hermitian integral ====
  recursive function coef_R1(zp,wpc,n,j) result(res)
    FIELD :: zp, wpc(3)
    integer, intent(in) :: n(3), j
    FIELD :: res, d2
    integer :: id_i(3), i
    FIELD :: fjs(0:j)

    res = 0
    d2 = dot_product(wpc, wpc)
        
    if(all(n==0)) then
       call igamma(j, zp*d2, fjs); check_err()
       res = (-2*zp)**j * fjs(j)
    else
       do i = 1, 3
          id_i(:) = 0; id_i(i) = 1
          if(n(i)>0) then
             res = wpc(i) * coef_R1(zp,wpc,n(:)-id_i(:),j+1)
             if(n(i)>1) then
                res = res + (n(i)-1) * coef_R1(zp,wpc,n(:)-2*id_i(:),j+1)
             end if
          end if
       end do
    end if
    
  end function coef_R1  
  subroutine coef_R_fast(zp,wpc,maxn, cr)
    FIELD, intent(in) :: zp
    FIELD, intent(in) :: wpc(3)
    integer, intent(in) :: maxn
    FIELD, intent(out) :: cr(0:,0:,0:)
    FIELD, allocatable :: Fj(:)
    FIELD, allocatable :: rmap(:,:,:,:)
    FIELD :: d2, tmp
    integer j, nnn, nx, ny, nz
    
    cr(:,:,:) = 0
    d2 = dot_product(wpc, wpc)

    allocate(Fj(0:maxn*3))
    allocate(rmap(0:maxn,0:maxn,0:maxn,0:3*maxn))

    call igamma(3*maxn, zp*d2, Fj(:)); check_err()

    tmp = 1
    do j = 0, 3*maxn
       rmap(0,0,0,j) = tmp * Fj(j)
       tmp = tmp*(-2*zp)
    end do

    do nnn = 1, 3*maxn
       do nx = 0, min(maxn, nnn)
          do ny = 0, min(maxn, nnn-nx)
             do nz = 0, min(maxn, nnn-nx-ny)
                do j = 0, 3*maxn-nnn
                   if(nx>0) then
                      tmp = wpc(1)*rmap(nx-1,ny,nz,j+1)
                      if(nx>1) then
                         tmp = tmp + (nx-1)*rmap(nx-2,ny,nz,j+1)
                      end if
                      rmap(nx,ny,nz,j) = tmp
                   else if(ny>0) then
                      tmp = wpc(2)*rmap(nx,ny-1,nz,j+1)
                      if(ny>1) then
                         tmp = tmp + (ny-1)*rmap(nx,ny-2,nz,j+1)
                      end if
                      rmap(nx,ny,nz,j) = tmp                      
                   else if(nz>0) then
                      tmp = wpc(3)*rmap(nx,ny,nz-1,j+1)
                      if(nz>1) then
                         tmp = tmp + (nz-1)*rmap(nx,ny,nz-2,j+1)
                      end if
                      rmap(nx,ny,nz,j) = tmp                      
                   end if
                end do
             end do
          end do
       end do
    end do

    cr(0:maxn,0:maxn,0:maxn) = rmap(:,:,:,0)
        
  end subroutine coef_R_fast
  subroutine coef_R(zp,wpc,maxn, cr)
    use Mod_Timer
    FIELD, intent(in) :: zp
    FIELD, intent(in) :: wpc(3)
    integer, intent(in) :: maxn
    FIELD, intent(out) :: cr(0:,0:,0:)
    integer nx, ny, nz, n(3)    

    select case(method_coef_R_)
    case(0)
       do nx = 0, maxn
          do ny = 0, maxn
             do nz = 0, maxn
                n(1) = nx; n(2) = ny; n(3) = nz
                cr(nx,ny,nz) = coef_R1(zp,wpc,n,0)
                check_err()
             end do
          end do
       end do
    case(1)
       call coef_R_fast(zp,wpc,maxn, cr); check_err()
    end select
    
  end subroutine coef_R
end module Mod_Molfunc
