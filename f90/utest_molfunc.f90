#include "macros.fpp"
#include "macros_utest.fpp"

! from enint/py/for_f90/igamma.py
!0 (0.1-0.2j) (0.963925032536+0.0626298665132j)  
!1 (0.1-0.2j) (0.311391276351+0.0371146636711j)
!2 (0.1-0.2j) (0.184211823869+0.026326381207j)
!3 (0.1-0.2j) (0.130521284118+0.0203848768131j)
!0 (0.3+0.1j) (0.907584447353-0.0279102336187j)
!1 (0.3+0.1j) (0.278724938324-0.0161613453311j)
!2 (0.3+0.1j) (0.161323588566-0.0113172327027j)
!3 (0.3+0.1j) (0.112937216594-0.00869198791243j)
!0 (21.2+15.5j) (0.164390140077-0.0536861042201j)
!1 (21.2+15.5j) (0.00192328174186-0.00267235467192j)
!2 (21.2+15.5j) (-1.40915382131e-06-0.000188051417904j)
!3 (21.2+15.5j) (-1.0673882707e-05-1.43718537129e-05j)
!0 (23+15.5j) (0.160935386952-0.0491668665616j)
!1 (23+15.5j) (0.00191058009153-0.00235640976908j)
!2 (23+15.5j) (1.44666718823e-05-0.0001634281764j)
!3 (23+15.5j) (-7.15113017015e-06-1.29446918596e-05j)

module Mod_TestMolfnc
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_Molfunc
  implicit none
contains
  subroutine TestMolfunc_run()

    call Utest_sub_begin("coef_d")
    call test_molfunc()
    call Utest_sub_end()
    
  end subroutine TestMolfunc_run
  subroutine test_molfunc()
    complex(kind(0d0)) :: ref(4*4), calc(0:3)
    complex(kind(0d0)) :: zs(4)
    integer iz, m, i
    double precision, parameter :: eps = 1.0d-10

    zs(:) = (/(0.1d0, -0.2d0), (0.3d0, 0.1d0), (21.2d0, 15.5d0), (23.0d0, 15.5d0)/)
    ref(:) = (/ &
    (0.963925032536d0, +0.0626298665132d0), &
    (0.311391276351d0, +0.0371146636711d0), &      
    (0.184211823869d0, +0.026326381207d0) , &      
    (0.130521284118d0, +0.0203848768131d0) ,&     
    (0.907584447353d0, -0.0279102336187d0), &  
    (0.278724938324d0, -0.0161613453311d0), &   
    (0.161323588566d0, -0.0113172327027d0), &  
    (0.112937216594d0, -0.00869198791243d0),& 
    (0.164390140077d0, -0.0536861042201d0), &
    (0.00192328174186d0, -0.00267235467192d0), &
    (-1.40915382131d-06, -0.000188051417904d0) ,&
    (-1.0673882707d-05, -1.43718537129d-05)  ,&
    (0.160935386952d0, -0.0491668665616d0)        ,&     
    (0.00191058009153d0, -0.00235640976908d0),&
    (1.44666718823d-05, -0.0001634281764d0) ,       &
    (-7.15113017015d-06, -1.29446918596d-05) /)

    i = 0
    do iz = 1, 4
       call igamma(3, zs(iz), calc); check_err()
       do m = 0, 3
          i = i+1
          call expect_near(ref(i), calc(m), eps)
       end do
    end do

  end subroutine test_molfunc
end module Mod_TestMolfnc
  
program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestMolfnc

  call utest_new
  call ErrHandle_new
  
  call TestMolfunc_run()

  call ErrHandle_delete
  call Timer_result
  call utest_delete
  
end program main
