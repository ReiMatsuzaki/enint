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

    call Utest_sub_begin("igamma_c")
    call test_igamma_c()
    call Utest_sub_end()

    call Utest_sub_begin("coef_d")
    call test_coef_d()
    call Utest_sub_end()

    call Utest_sub_begin("igamma_r")
    call test_igamma_r()
    call Utest_sub_end()

    call Utest_sub_begin("igamma_rf")
    call test_igamma_r_fast()
    call Utest_sub_end()    

    call Utest_sub_begin("coef_r")
    call test_coef_r()
    call Utest_sub_end(); check_err()

    call Utest_sub_begin("fast_cr")
    call test_coef_r_fast()
    call Utest_sub_end(); check_err()
    
  end subroutine TestMolfunc_run
  subroutine test_igamma_c()
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
       call igamma_c(3, zs(iz), calc); check_err()
       do m = 0, 3
          i = i+1
          call expect_near(ref(i), calc(m), eps)
          if(get_err().ne.0) then
             write(0,*) "z = ", zs(iz)
             write(0,*) "m = ", m
             return
             !call ErrHandle_ierr0
          end if
       end do
    end do

  end subroutine test_igamma_c
  subroutine test_coef_d()
    
    ![[[ 1.      0.      0.      0.      0.      0.    ]
    ![-0.2     0.5     0.      0.      0.      0.    ]
    ![ 0.54   -0.2     0.25    0.      0.      0.    ]]
    !
    ![[-0.1     0.5     0.      0.      0.      0.    ]
    ![ 0.52   -0.15    0.25    0.      0.      0.    ]
    ![-0.254   0.79   -0.125   0.125   0.      0.    ]]
    !
    ![[ 0.51   -0.1     0.25    0.      0.      0.    ]
    ![-0.202   0.775  -0.1     0.125   0.      0.    ]
    ![ 0.8154 -0.456   0.7825 -0.075   0.0625  0.    ]]]
    
    call expect_eq(1.0d0,  coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,0,0))
    call expect_eq(-0.2d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,1,0))
    call expect_eq(0.5d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,1,1))
    call expect_eq(0.54d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,2,0))
    
  end subroutine test_coef_d
  subroutine test_igamma_r()
    ! see ${ENINT}/py/for_f90/
    !    0 1.0 0.746824132812
    !    1 1.1 0.179771018184
    !    2 3.3 0.0251274571212
    !    3 0.0 0.142857142857
    double precision :: fs(0:5)
    
    call igamma_r(2, 1.0d0, fs); check_err()
    call expect_eq(0.746824132812d0, fs(0))

    call igamma(2, 1.1d0, fs); check_err()
    call expect_eq(0.179771018184d0, fs(1))

    call igamma_r(2, 3.3d0, fs); check_err()
    call expect_eq(0.0251274571212d0, fs(2))

    call igamma_r(3, 0.0d0, fs); check_err()
    call expect_eq(0.142857142857d0, fs(3))
    
  end subroutine test_igamma_r
  subroutine test_igamma_r_fast()
    ! see ${ENINT}/py/for_f90/
    !    0 1.0 0.746824132812
    !    1 1.1 0.179771018184
    !    2 3.3 0.0251274571212
    !    3 0.0 0.142857142857
    double precision :: f0(0:10), f1(0:10)
    integer m
    integer, parameter :: nz = 4
    double precision :: z(nz) = (/1.0d0, 0.0d0, 0.00001d0, 13.0d0/)
    integer i
    do i = 1, nz
       method_igamma_ = 0
       call igamma(8, z(i), f0); check_err()
       method_igamma_ = 1
       call igamma(8, z(i), f1); check_err()
       do m = 0, 6
          call expect_not_NaN_d(f1(m))
          if(get_err().eq.1) then
             begin_err(1)
             write(*,*) "m=", m
             write(*,*) "z=", z(i)
             end_err()
          end if
          call assert_near_d(f0(m), f1(m), 1.0d-13)
          if(get_err().eq.1) then
             begin_err(1)
             write(*,*) "m=", m
             write(*,*) "z=", z(i)
             end_err()
          end if
       end do
    end do
    
    
  end subroutine test_igamma_r_fast
  subroutine test_coef_r()
    ! see ${ENINT}/py/for_f90
    ! print coef_R_list(1.1, [0.0,0.1,0.2], [0.2,0.3,0.4], 1, 0)    
    ![[[ 0.95768901  0.13558003]
    !  [ 0.13558003  0.0352501 ]]
    !
    ! [[ 0.13558003  0.0352501 ]
    !  [ 0.0352501   0.0109848 ]]]
    double precision :: cr(0:1,0:1,0:1)
    double precision :: wp(3), wc(3), wpc(3)

    wp(:) = (/0.0d0,0.1d0,0.2d0/)
    wc(:) = (/0.2d0,0.3d0,0.4d0/)
    wpc(:) = wp(:)-wc(:)
    call coef_R(1.1d0, wpc(:), 1, cr)

    call expect_eq(0.95768901d0, cr(0,0,0))
    call expect_eq(0.0352501d0, cr(1,1,0))
    call expect_eq(0.0352501d0, cr(0,1,1))
    
  end subroutine test_coef_r
  subroutine test_coef_r_fast()
    double precision :: cr0(0:3,0:3,0:3), cr1(0:3,0:3,0:3)
    double precision :: zp, wp(3), wc(3), wpc(3)
    integer :: x,y,z

    zp = 1.1d0
    wp(:) = (/0.0d0,0.1d0,0.2d0/)
    wc(:) = (/0.3d0,0.0d0,0.1d0/)
    wpc(:) = wp(:)-wc(:)
    method_coef_R_ = 0
    call coef_R(zp, wpc, 3, cr0); check_err()
    method_coef_R_ = 1
    call coef_R(zp, wpc, 3, cr1); check_err()

    do x = 0, 3
       do y = 0, 3
          do z = 0, 3
             call assert_near_d(cr0(x,y,z), cr1(x,y,z), 1.0d-7)
             if(get_err().ne.0) then
                begin_err(1)
                write(0,*) "x,y,z=", x,y,z
                end_err()
                get_err() = 0
             end if
          end do
       end do
    end do
    
  end subroutine test_coef_r_fast
  
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
