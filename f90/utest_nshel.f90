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
