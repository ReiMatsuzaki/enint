#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_TestNshel
  use Mod_ErrHandle
  use Mod_UtestCheck
  use Mod_Nshel
  implicit none
contains
  subroutine TestNshel_run()
    call test_coef()
    call test_smat()
    call test_h2()
  end subroutine TestNshel_run
  subroutine test_coef()
    write(*,*) "--------------------"
    write(*,*) "TestNshel_coef begin"
    
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
    call expect_eq(0.54d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,2,0))

    write(*,*) "TestNshel_coef end"
    write(*,*) "--------------------"
    
  end subroutine test_coef
  subroutine test_smat()
    type(Obj_Nshel) nshel
    integer :: ns(1,3)
    double precision :: coef_l(0:3,20)
    double precision, allocatable :: mat(:,:)

    write(*,*) "--------------------"
    write(*,*) "TestNshel_smat begin"

    call Nshel_new(nshel, 2, 3); check_err()
    nshel%nucs%ws(1,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%ws(2,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%zs(:) = (/1.0d0, 2.0d0/)
    ns = 0
    coef_l = 0
    coef_l(0,1) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/), 1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"s"/), 1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"s"/), 1, (/1.4d0/), coef_l, 2); check_err()
    call Nshel_setup(nshel); check_err()

    allocate(mat(nshel%nbasis, nshel%nbasis))
    call Nshel_s(nshel, mat); check_err()
    
    call Nshel_delete(nshel); check_err()

    write(*,*) "TestNshel_smat end"
    write(*,*) "--------------------"

  end subroutine test_smat
  subroutine test_h2()
    use Mod_math
    integer :: ifile = 12323
    double precision, allocatable :: calc(:,:), ref(:,:)
    type(Obj_Nshel) nshel
    integer :: num

    write(*,*) "--------------------"
    write(*,*) "TestNshel_h2 begin"

    ! -- new --
    call Nshel_new_file(nshel, "../gms/h2/out/nshel.json"); check_err()
    call Nshel_setup(nshel); check_err()
    call Nshel_dump(nshel); check_err()
    num = nshel%nbasis
    allocate(ref(num,num), calc(num,num))

    ! -- s matrix --
    call Nshel_s(nshel, calc); check_err()
    call open_r(ifile, "../gms/h2/out/s.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()
    call expect_near_dmat(ref, calc, 10.0d0**(-7))
    close(ifile)

    ! -- t matrix --
    call Nshel_t(nshel, calc); check_err()
    call open_r(ifile, "../gms/h2/out/t.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()

    call expect_near_dmat(ref, calc, 10.0d0**(-7))
    close(ifile)
    
    ! -- delete --
    call Nshel_delete(nshel); check_err()
    deallocate(calc, ref)
    close(ifile)

    write(*,*) "TestNshel_h2 end"
    write(*,*) "--------------------"    
    
  end subroutine test_h2
end module Mod_TestNshel

program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestNshel

  call utest_new
  call ErrHandle_new

  call TestNshel_run

  call ErrHandle_delete
  call utest_delete
  
end program main
