#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_TestNshel
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_Nshel
  implicit none
contains
  subroutine TestNshel_run()
    
    call Utest_sub_begin("coef_d")
    call test_coef_d()
    call Utest_sub_end()

    call Utest_sub_begin("gammainc")
    call test_gammainc()
    call Utest_sub_end()

    call Utest_sub_begin("gi_f")
    call test_gammainc_fast()
    call Utest_sub_end()

    call Utest_sub_begin("coef_r")
    call test_coef_r()
    call Utest_sub_end(); check_err()

    call Utest_sub_begin("fast_cr")
    call test_coef_r_fast()
    call Utest_sub_end(); check_err()

    call Utest_sub_begin("smat")
    call test_smat()
    call Utest_sub_end()

    call Utest_sub_begin("h2")
    call test_h2()
    call Utest_sub_end()

    call Utest_sub_begin("hcp")
    call test_hcp()
    call Utest_sub_end()
    
  end subroutine TestNshel_run
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
  subroutine test_gammainc()
    ! see ${ENINT}/py/for_f90/
    !    0 1.0 0.746824132812
    !    1 1.1 0.179771018184
    !    2 3.3 0.0251274571212
    !    3 0.0 0.142857142857
    double precision :: fs(0:5)
    
    call mole_gammainc(2, 1.0d0, fs); check_err()
    call expect_eq(0.746824132812d0, fs(0))

    call mole_gammainc(2, 1.1d0, fs); check_err()
    call expect_eq(0.179771018184d0, fs(1))

    call mole_gammainc(2, 3.3d0, fs); check_err()
    call expect_eq(0.0251274571212d0, fs(2))

    call mole_gammainc(3, 0.0d0, fs); check_err()
    call expect_eq(0.142857142857d0, fs(3))
    
  end subroutine test_gammainc
  subroutine test_gammainc_fast()
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
       call mole_gammainc(8, z(i), f0, 0); check_err()
       call mole_gammainc(8, z(i), f1, 1); check_err()
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
    
    
  end subroutine test_gammainc_fast
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
    call coef_R(zp, wpc, 3, cr0, 0); check_err()
    call coef_R(zp, wpc, 3, cr1, 1); check_err()

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
  subroutine test_smat()
    type(Obj_Nshel) nshel
    integer :: ns(1,3)
    double precision :: coef_l(0:3,20)
    double precision, allocatable :: mat(:,:)

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

  end subroutine test_smat
  subroutine test_h2()    
    use Mod_math
    integer :: ifile = 12323
    double precision, allocatable :: calc(:,:), ref(:,:), v(:,:)
    type(Obj_Nshel) nshel
    integer :: num

    ! -- new --
    call Nshel_new_file(nshel, "../gms/h2/out/nshel.json"); check_err()
    call Nshel_setup(nshel); check_err()
    num = nshel%nbasis
    allocate(ref(num,num), calc(num,num), v(num,num))

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

    ! -- h matrix --
    call Nshel_v(nshel, v); check_err()
    call Nshel_t(nshel, calc); check_err()
    calc(:,:) = calc(:,:) + v(:,:)
    call open_r(ifile, "../gms/h2/out/h.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()    
    call expect_near_dmat(ref, calc, 10.0d0**(-7))
    
    ! -- delete --
    call Nshel_delete(nshel); check_err()
    deallocate(calc, ref)
    close(ifile)
    
  end subroutine test_h2
  subroutine test_hcp()
    use Mod_math
    integer :: ifile = 12323
    double precision, allocatable :: calc(:,:), ref(:,:), v(:,:)
    type(Obj_Nshel) nshel
    integer :: num

    ! -- new --
    call Nshel_new_file(nshel, "../gms/hcp/out/nshel.json"); check_err()
    call Nshel_setup(nshel); check_err()
    num = nshel%nbasis
    allocate(ref(num,num), calc(num,num), v(num,num))

    ! -- s matrix --
    call Nshel_s(nshel, calc); check_err()
    call open_r(ifile, "../gms/hcp/out/s.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()
    call expect_prop_dmat(calc, "overlap"); check_err()
    call expect_near_dmat(ref, calc, 10.0d0**(-9)); check_err()
    close(ifile)

    ! -- t matrix --
    call Nshel_t(nshel, calc); check_err()
    call open_r(ifile, "../gms/hcp/out/t.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()
    call expect_near_dmat(ref, calc, 10.0d0**(-9))

    ! -- h matrix --
    call Nshel_v(nshel, v); check_err()
    call Nshel_t(nshel, calc); check_err()
    calc(:,:) = calc(:,:) + v(:,:)
    call open_r(ifile, "../gms/hcp/out/h.csv"); check_err()    
    call load_dmat(ifile, ref); check_err()
    call expect_prop_dmat(calc, "hermite"); check_err()
    call expect_near_dmat(ref, calc, 10.0d0**(-9)); check_err()
    
    ! -- delete --
    call Nshel_delete(nshel); check_err()
    deallocate(calc, ref)
    close(ifile)
    
  end subroutine test_hcp
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
  call Timer_result
  call utest_delete
  
end program main
