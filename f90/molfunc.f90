#include "macros.fpp"

module Mod_Molfunc
  use Mod_ErrHandle
  implicit none
contains
  subroutine igamma_f1(maxm, z, res)
    integer, intent(in) :: maxm
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(:)
    double precision, parameter :: esp = 1.0d-10
    integer, parameter :: NF_MAX=50, NF_0=10
    double precision :: anR(0:maxm), anI(0:maxm), z2, az
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

    anR(0) = -x/(2*z2)*exp(-x)*cos(y) + y/(2*z2)*exp(-x)*sin(y)
    anI(0) = +y/(2*z2)*exp(-x)*cos(y) + x/(2*z2)*exp(-x)*sin(y)
    do n = 1, NF_MAX
       anR(n) = -(2*n-1) * (x/(2*z2)*anR(n-1) + y/(2*z2)*anI(n-1))
       anI(n) = +(2*n-1) * (y/(2*z2)*anR(n-1) - x/(2*z2)*anI(n-1))
       if(n>NF_0) then
          if(c*delta > abs(anR(n)) and s*delta>abs(anI(n))) then
             convq = True
             exit
          end if
       end if
    end do
    if(.not. convq) then
       throw_err("not converged",1)
    end if
    
    fmR[0] = +c + sum(anR)
    fmI[0] = -s + sum(anI)
    do m = 1, maxm
       fmR[m] = (2*m-1) * (x/(2*z2)*fmR(m-1) + y/(2*z2)*fmI(m-1))
       fmI[m] = (2*m-1) * (x/(2*z2)*fmI(m-1) - y/(2*z2)*fmR(m-1))
    end do

    bmR(0) = 0
    bmI(0) = 0
    do m = 1, maxm
       bmR(m) = anR(0) + (2*m-1)*(x/(2*z2)*bmR(m-1) + y/(2*z2)*bmI(m-1))
       bmI(m) = anI(0) + (2*m-1)*(x/(2*z2)*bmI(m-1) - y/(2*z2)*bmR(m-1))
    end do

    res(:) = (fmR(:) + bmR(:)) + ii * (fmI(:)+bmI(:))
    
  end subroutine igamma_f1
  subroutine igamma_f2(maxm, z, res)
    integer, intent(in) :: maxm
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(:)
    double precision, parameter :: esp = 1.0d-10
    integer, parameter NR=47
    complex(kind(0d0)) :: bn(0:NR), an(0:NR)
    double precision tmp1, tmp3, tmp5, tmp7

    x = real(z)
    y = aimag(z)

    if(x < -eps) then
       throw_err("Re[z] must be positive", 1)
    end if
    if(y < -eps) then
       call igamma_c(maxm, z, res); check_err()
       res(:) = conjg(res(:))
       return 
    end if

    bn(0) = 1
    bn(1) = 1 + z/2
    do n = 2, NR
       bn[n] = bn[n-1] + z*z/(4*(2*n-1)*(2*n-3))*bn[n-2]
    end do

    do m = 0, maxm
       an(0) = 1
       tmp1 = 2*m+1; tmp3 = tmp1+2; tmp5 = tmp3+2; tmp7 = tmp5+2       
       t1 = tmp1/tmp3
       t2 = tmp1/(tmp3*tmp5)
       t3 = tmp1*(tmp1**2+44)/(60*tmp3*tmp5*tmp7)
       an[1] = bn[1] -t1*z
       an[2] = bn[2] -t1*z -t2*z**2
       an[3] = bn[3] -t1*z -t2*z**2  -t3*z**3
       do n = 4, NR
          f1 = (2.0d0*n-2*m-5)/(2*(2*n-3)*(2*n+2*m+1))
          f2 = 1.0d0/(4*(2*n-1)*(2*n-3))
          f3 = -f1/(4*(2*n-3)*(2*n-5))
          e = -f1
          an[n] = (1+f1*z)*an(n-1) + (e+f2*z)*z*an(n-2) + f3*z**3*an(n-3)
       end do
       res(m) = 1.0d0/(2*m+1) * an(NR)/bn(NR)
    end do
    
  end subroutine igamma_f2
end module Mod_Molfunc
