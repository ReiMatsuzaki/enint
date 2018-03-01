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
      
!    call Utest_sub_begin("smat")
!    call test_smat()
!    call Utest_sub_end()
!
!    call Utest_sub_begin("ao_at")
!    call test_ao_at()
!    call Utest_sub_end()
!
!    call Utest_sub_begin("ao_at_d")
!    call test_ao_at_d()
    !    call Utest_sub_end()

!    call Utest_sub_begin("ao_at_2")
!    call test_ao_at_2
!    call Utest_sub_end()

!    call Utest_sub_begin("aoao_at_s")
!    call test_aoao_at_s
    !    call Utest_sub_end()

    call Utest_sub_begin("aoao_at")
    call test_aoao_at
    call Utest_sub_end()
    
!    call Utest_sub_begin("ao_dw")
!    call test_dw()
!    call Utest_sub_end()

 !   call Utest_sub_begin("h2")
 !   call test_h2()
 !   call Utest_sub_end()
!
!    call Utest_sub_begin("hcp")
!    call test_hcp()
!    call Utest_sub_end()
    
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
  subroutine test_ao_at()
    
    type(Obj_Nshel) nshel
    double precision :: coef_l(0:3,20)
    double precision :: rs(2,3)
    double precision, allocatable :: ao_rs(:,:), ao_rs2(:,:)
    integer :: num_basis
    double precision :: ref, d

    d = 1.2d0

    call Nshel_new(nshel, 2, 3); check_err()
    nshel%nucs%ws(1,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%ws(2,:) = d*(/1.0d0,0.0d0,0.0d0/);
    nshel%nucs%zs(:)   = (/1.0d0, 0.3d0/)
    coef_l = 0
    coef_l(0,1) = 1.0d0
    coef_l(1,1) = 1.0d0
    coef_l(2,1) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/),   1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"p"/),   1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"d"/), 1, (/1.2d0/), coef_l, 2); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    
    rs = reshape((/0.0d0, 1.3d0, &
         0.0d0, 0.1d0, &
         0.0d0, -0.1d0/), (/2,3/))
    call Nshel_num_basis(nshel, num_basis); check_err()
    allocate(ao_rs(num_basis, size(rs,1)))
    allocate(ao_rs2(num_basis, size(rs,1)))
    call Nshel_ao_at(nshel, rs, (/0,0,0/), ao_rs(:,:)); check_err()

    call expect_eq(nshel%shels(1)%coef(1,1), ao_rs(1,1)); check_err()
    
    ref = -nshel%shels(2)%coef(1,1)*d*exp(-1.3d0*d**2)
    call expect_eq(ref, ao_rs(2,1)); check_err()

    call expect_eq(0.0d0, ao_rs(3,1))
    call expect_eq(0.0d0, ao_rs(4,1))

    ref = nshel%shels(3)%coef(1,1) * d**2 * exp(-1.2d0*d*d)    
    call expect_eq(ref, ao_rs(5,1))
    
  end subroutine test_ao_at
  subroutine test_ao_at_2()
    
    type(Obj_Nshel) nshel
    double precision :: coef_l(0:3,20)
    double precision, parameter :: xmax = 5.0d0
    integer, parameter :: n1 = 50
    integer, parameter :: nr = n1**3
    double precision :: rs(nr,3)
    double precision, allocatable :: ao_rs(:,:)
    integer :: num_basis, i, j, k, idx
    double precision :: dx
    
    call Nshel_new(nshel, 2, 3); check_err()
    nshel%nucs%ws(1,:) = (/0.1d0,0.2d0,-0.1d0/);
    nshel%nucs%ws(2,:) = (/0.3d0,-0.1d0,0.2d0/);   ! location 
    nshel%nucs%zs(:)   = (/1.0d0, 0.3d0/)         ! charge
    coef_l = 0
    coef_l(0,1) = 1.0d0
    coef_l(1,1) = 1.0d0
    coef_l(2,1) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/),   1, (/1.0d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"p"/),   1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"d"/), 1, (/1.2d0/), coef_l, 2); check_err()
    call Nshel_setup(nshel, .true.); check_err()


    dx = 2*xmax/n1
    idx = 0
    do i = 1, n1
       do j = 1, n1
          do k = 1, n1
             idx = idx + 1
             rs(idx, 1) = (i - n1/2)*dx
             rs(idx, 2) = (j - n1/2)*dx
             rs(idx, 3) = (k - n1/2)*dx
          end do
       end do
    end do
    
    call Nshel_num_basis(nshel, num_basis); check_err()
    allocate(ao_rs(num_basis, nr))
    call Nshel_ao_at(nshel, rs, (/0,0,0/), ao_rs(:,:)); check_err()

    call expect_near(1.0d0, sum(ao_rs(1, :)**2 * dx**3), 1.0d10)
    
    call expect_near(1.0d0, sum(ao_rs(2, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(3, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(4, :)**2) * dx**3, 1.0d-10)
    
    call expect_near(1.0d0, sum(ao_rs(5, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(6, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(7, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(8, :)**2) * dx**3, 1.0d-10)
    call expect_near(1.0d0, sum(ao_rs(9, :)**2) * dx**3, 1.0d-10)
    
  end subroutine test_ao_at_2
  subroutine test_ao_at_d()
    type(Obj_Nshel) nshel
    double precision :: coef_l(0:3,20)
    integer :: num_basis, mu
    double precision :: rs(7,3), r0(3), dx, ref
    double precision, allocatable :: phi(:,:), phi_x(:,:), phi_y(:,:), phi_z(:,:)

    dx = 0.001d0

    call Nshel_new(nshel, 1, 3); check_err()
    !    nshel%nucs%ws(1,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%ws(1,:) = (/0.1d0,0.2d0,0.3d0/);
    coef_l = 0
    coef_l(0,1) = 1.0d0
    coef_l(1,1) = 1.0d0
    coef_l(2,1) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/), 1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"p"/), 1, (/1.3d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 3, (/"d"/), 1, (/1.2d0/), coef_l, 1); check_err()
    call Nshel_setup(nshel, .true.); check_err()

    r0 = (/0.3d0, 0.2d0, 0.1d0/)
    rs(1,:) = r0
    rs(2,:) = r0 + dx * (/1.0d0, 0.0d0, 0.0d0/)
    rs(3,:) = r0 - dx * (/1.0d0, 0.0d0, 0.0d0/)
    rs(4,:) = r0 + dx * (/0.0d0, 1.0d0, 0.0d0/)
    rs(5,:) = r0 - dx * (/0.0d0, 1.0d0, 0.0d0/)
    rs(6,:) = r0 + dx * (/0.0d0, 0.0d0, 1.0d0/)
    rs(7,:) = r0 - dx * (/0.0d0, 0.0d0, 1.0d0/)

    call Nshel_num_basis(nshel, num_basis)
    allocate(phi(  num_basis, size(rs, 1)))
    allocate(phi_x(num_basis, size(rs, 1)))
    allocate(phi_y(num_basis, size(rs, 1)))
    allocate(phi_z(num_basis, size(rs, 1)))
    call Nshel_ao_at(nshel, rs(:,:), (/0,0,0/), phi(  :,:)); check_err()
    call Nshel_ao_at(nshel, rs(:,:), (/1,0,0/), phi_x(:,:)); check_err()
    call Nshel_ao_at(nshel, rs(:,:), (/0,1,0/), phi_y(:,:)); check_err()
    call Nshel_ao_at(nshel, rs(:,:), (/0,0,1/), phi_z(:,:)); check_err()

    do mu = 1, num_basis
       ref = (phi(mu,2)-phi(mu,3))/(2*dx)
       call expect_near(ref, phi_x(mu, 1), 1.0d-5); check_err()
       ref = (phi(mu,4)-phi(mu,5))/(2*dx)
       call expect_near(ref, phi_y(mu, 1), 1.0d-5); check_err()
       ref = (phi(mu,6)-phi(mu,7))/(2*dx)
       call expect_near(ref, phi_z(mu, 1), 1.0d-5)
       if(get_err() .ne. 0) then
          write(0,*) "mu:", mu
          return
       end if       
    end do
        
  end subroutine test_ao_at_d
  subroutine test_aoao_at_s()
    use Mod_const, only : pi
    type(Obj_Nshel)  :: nshel
    double precision :: coef_l(0:3,20)
    integer, parameter :: natom=2, num_sh=2
    integer :: nao
    double precision :: zs(1), zeta1, zeta2
    double precision, allocatable :: aoao(:,:,:)
    double precision :: ref, d1, d2

    call Nshel_new(nshel, natom, num_sh)
    nshel%nucs%ws(1,:) = (/0.0d0, 0.0d0, 0.1d0/)
    nshel%nucs%ws(2,:) = (/0.0d0, 0.0d0, 0.2d0/)
    coef_l = 0
    coef_l(0,1) = 1
    coef_l(1,1) = 1
    coef_l(2,1) = 1
    call Nshel_set(nshel, 1, (/"s"/), 1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"d"/), 1, (/1.1d0/), coef_l, 1); check_err()    
    call Nshel_setup(nshel, .true.); check_err()
    !call Nshel_dump(nshel)
    !   j  jj   n      
    !   2   1  2  0  0     1.94467
    !   3   2  0  2  0     1.94467
    !   4   3  0  0  2     1.94467
    !   5   4  1  1  0     3.36827
    !   6   5  1  0  1     3.36827
    !   7   6  0  1  1     3.36827

    call Nshel_num_basis(nshel, nao)
    allocate(aoao(nao, nao, 1))
    
    zs(1) = 0.3d0
    call Nshel_aoao_at_intxy(nshel, zs, 0, 0, aoao); check_err()

    zeta1 = nshel%shels(1)%zeta(1)
    zeta2 = nshel%shels(2)%zeta(1)
    d1    = nshel%shels(1)%w(3)
    d2    = nshel%shels(2)%w(3)
    
    call expect_near(0.0d0, aoao(1,5,1), 1.0d-8); check_err()
    call expect_near(0.0d0, aoao(1,6,1), 1.0d-8); check_err()
    call expect_near(0.0d0, aoao(1,7,1), 1.0d-8); check_err()

    ref = sqrt(2.0d0*zeta1/pi) * exp(-2*zeta1*(zs(1)-d1)**2)
    call expect_near(ref, aoao(1,1,1), 1.0d-8); check_err()
    
    ref = sqrt(2.0d0*zeta2/pi) * exp(-2*zeta2*(zs(1)-d2)**2)
    call expect_near(ref, aoao(2,2,1), 1.0d-8); check_err()
    call expect_near(ref, aoao(3,3,1), 1.0d-8); check_err()    
    
  end subroutine test_aoao_at_s
  subroutine test_aoao_at()
    type(Obj_Nshel)  :: nshel
    double precision :: coef_l(0:3,20)
    integer, parameter :: natom=2, num_sh=3
    double precision, parameter :: maxx=4.5d0
    integer, parameter :: nx=1500, ny=nx, nz=1
    double precision, parameter :: dr=maxx/nx
    integer :: ix,iy, iz, nao, idx
    double precision :: zs(nz), rs(nx*ny*nz,3)
    
    double precision, allocatable :: ao(:,:), aoao(:,:,:)
    double precision :: ref

    call Nshel_new(nshel, natom, num_sh)
    nshel%nucs%ws(1,:) = (/0.0d0, 0.0d0, 0.0d0/)
    nshel%nucs%ws(2,:) = (/0.0d0, 0.0d0, 1.2d0/)
    nshel%nucs%zs(:)   = (/1.0d0, 3.0d0/)
    coef_l = 0
    coef_l(0,1) = 1
    coef_l(1,1) = 1
    coef_l(2,1) = 1
    call Nshel_set(nshel, 1, (/"s"/), 1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"p"/), 1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"d"/), 1, (/1.2d0/), coef_l, 2); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    
    do iz = 1, nz
       zs(iz) = 0.3d0
       do ix = 1, nx
          do iy = 1, ny
             idx = (iz-1)*nx*ny + (iy-1)*nx + ix
             rs(idx,1) = dr*(ix-1)
             rs(idx,2) = dr*(iy-1)
             rs(idx,3) = zs(iz)
          end do
       end do
    end do

    call Nshel_num_basis(nshel, nao)
    allocate(ao(nao, nx*ny*nz))
    call Nshel_ao_at(nshel, rs, (/0,0,0/), ao); check_err()
    allocate(aoao(nao, nao, nz))
    call Nshel_aoao_at_intxy(nshel, zs, 0, 0, aoao); check_err()

    ref = sum(ao(1,:)*ao(1,:))*dr*dr*4
    call expect_near(ref, aoao(1,1,1), 0.004d0); check_err()

    ref = sum(ao(5,:)*ao(5,:))*dr*dr*4
    call expect_near(ref, aoao(5,5,1), 0.001d0); check_err()

    ref = sum(ao(1,:)*ao(5,:))*dr*dr*4
    call expect_near(ref, aoao(1,5,1), 0.001d0); check_err()
    
  end subroutine test_aoao_at
  subroutine test_dw()
    type(Obj_Nshel) nshel
    integer :: ns(1,3)
    double precision :: coef_l(0:3,20)
    integer :: num
    double precision, allocatable :: M(:,:,:)
    integer :: i, j
    
    call Nshel_new(nshel, 2, 3); check_err()
    nshel%nucs%ws(1,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%ws(2,:) = (/0.2d0,0.1d0,0.3d0/);
    nshel%nucs%zs(:) = (/1.0d0, 2.0d0/)
    ns = 0
    coef_l = 0
    coef_l(0,1:2) = (/1.0d0, 0.5d0/)
    coef_l(1,1) = 1.0d0
    coef_l(2,1:3) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/), 2, (/1.1d0,1.0d0/), coef_l, 1)
    check_err()
    call Nshel_set(nshel, 2, (/"p"/), 1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"d"/), 1, (/1.3d0/), coef_l, 2); check_err()
    !call Nshel_set(nshel, 3, (/"p"/), 1, (/1.4d0/), coef_l, 2); check_err()
    
    !call Nshel_new_file(nshel, "../gms/hcp/out/nshel.json"); check_err()
    call Nshel_setup(nshel); check_err()
    !call Nshel_dump(nshel)
    
    call Nshel_num_basis(nshel, num)    
    allocate(M(3,num,num))
    call Nshel_dw(nshel, M)
    do i = 1, num
       do j = 1, num
          if(abs(M(1,i,j)+M(1,j,i)) > 1.0d-5) then
             write(*,*) i,j,M(1,i,j)
          end if
       end do
    end do
    write(*,*) 
    call expect_near(0.0d0, sum(M(1,:,:)+transpose(M(1,:,:))), 1.0d-8); check_err()
    
  end subroutine test_dw
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
