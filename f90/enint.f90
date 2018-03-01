#include "macros.fpp"

module Mod_Enint  
  use Mod_ErrHandle
  use Mod_Nshel
  use Mod_Fjson
  use Mod_Timer
  implicit none  
contains
  subroutine Enint_run
    character(100) :: runtype

    call ErrHandle_new()
    call Timer_new("Enint", .false.); check_err()

    write(*,*) " ==== Enint ===="

    if(iargc()<1) then
       begin_err(1)
       call print_help
       end_err()
    end if

    call getarg(1, runtype)
    write(*,*) "runtype:", trim(runtype)

    select case(runtype)
    case("na_rs")
       call na_rs_run; check_err()
    case("grid")
       call Grid_run; check_err()
    case("ao_grid")
       call AOGrid_run; check_err()
    case("aoao_grid")
       call AOAOGrid_run; check_err()       
    case("ao_mat")
       call AOMat_run; check_err()
    case default
       throw_err("unsupported",1)
    end select
  end subroutine Enint_run
  ! -- na_rs : NA matrix for nuclear=DVR --
  ! compute electron nuclear interaction matrix when
  ! electron part is represented by GTO and
  ! nuclear part is represented by DVR
  subroutine na_rs_run
    double precision, allocatable :: v(:,:)
    integer :: num, a, i, j, iargc
    type(Obj_Nshel) :: nshel
    integer, parameter :: iin=123, iout=124
    character(100) :: fn

    if(iargc()<4) then
       begin_err(1)
       call print_help
       end_err()
    end if
    
    call getarg(2, fn)
    write(*,*) "nshel json file:", fn
    call Nshel_new_file(nshel, fn); check_err()
    call Nucs_delete(nshel%nucs); check_err()
    call nucs_new(nshel%nucs, 1); check_err()    
    nshel%nucs%zs(1) = 1.0d0
    call Nshel_setup(nshel); check_err()
    num = nshel%nbasis
    allocate(v(num, num))

    call getarg(3, fn)
    write(*,*) "input rs file:", fn
    call open_r(iin, fn); check_err()

    call getarg(4, fn)
    write(*,*) "output file: ", fn
    call open_w(iout, fn); check_err()

    ! - start read -
    read(iin, *)
    write(iout,*) "a,i,j,val"
    a = 0
    do
       a = a+1
       
       read(iin, *, end=100) nshel%nucs%ws(1,:)
       
       call Nshel_v(nshel, v(:,:)); check_err()

       do i = 1, num
          do j = 1, num
             write(iout,'(i0,",",i0,",",i0,f20.10)') a, i, j, v(i,j)
          end do
       end do
    end do

100 continue
    
    close(iin)
    close(iout)
    call Nshel_delete(nshel); check_err()
    call Timer_result()
    call Timer_delete()
    call ErrHandle_delete()
    
  end subroutine na_rs_run
  ! -- Grid --
  subroutine Grid_run
    use Mod_ArgParser
    use Mod_StrUtil
    character(100) :: line, fn_out
    integer ix, iy, iz, nx, ny, nz
    double precision, allocatable :: xs(:), ys(:), zs(:)
    integer, parameter :: ifile = 411
    
    call arg_parse_s("-x", line); check_err()
    call str2vec(line, nx, xs, .false.)
    if(get_err().ne.0) then
       throw_err("error on converting -x", 1)
    end if
    allocate(xs(nx))
    call str2vec(line, nx, xs, .true.); check_err()

    call arg_parse_s("-y", line); check_err()
    call str2vec(line, ny, ys, .false.); check_err()
    if(get_err().ne.0) then
       throw_err("error on converting -y", 1)
    end if
    allocate(ys(ny))
    call str2vec(line, ny, ys, .true.); check_err()

    call arg_parse_s("-z", line); check_err()
    call str2vec(line, nz, zs, .false.); check_err()
    if(get_err().ne.0) then
       throw_err("error on converting -z", 1)
    end if
    allocate(zs(nz))
    call str2vec(line, nz, zs, .true.); check_err()

    call arg_parse_s("-o", fn_out); check_err()

    call open_w(ifile, fn_out); check_err()
    write(ifile, '(A)') "x,y,z"
    do ix = 1, nx
       do iy = 1, ny
          do iz = 1, nz
             write(ifile, '(f20.10,",",f20.10,",",f20.10)') xs(ix),ys(iy),zs(iz)
          end do
       end do
    end do
    close(ifile)
    
  end subroutine Grid_run
  subroutine AOGrid_run
    use Mod_fjson
    use Mod_ArgParser
    type(Obj_Nshel) :: nshel
    character(100) :: fn_nshel, fn_in, fn_out, str_type
    double precision :: rs(1,3)
    integer :: ir, num_basis, mu
    integer :: n_dr(3)
    double precision, allocatable :: ao(:,:)
    integer, parameter :: ifile_in = 192
    integer, parameter :: ifile_out = 193

    write(*,*) "AOGrid_run"
    
    call arg_parse_s("-t", str_type); check_err()
    call arg_parse_s("-i", fn_in); check_err()
    call arg_parse_s("-o", fn_out); check_err()
    call arg_parse_s("-nshel", fn_nshel); check_err()

    select case(str_type)
    case("000")
       n_dr = (/0,0,0/)
    case("100")
       n_dr = (/1,0,0/)
    case("010")
       n_dr = (/0,1,0/)
    case("001")
       n_dr = (/0,0,1/)
    case default
       throw_err("unsupported str_type", 1)
    end select

    write(*,*) "n_dr:", n_dr

    call open_r(ifile_in, fn_in); check_err()
    call open_w(ifile_out, fn_out); check_err()
    
    write(*,*) "nshel json file:", trim(fn_nshel)
    call Nshel_new_file(nshel, fn_nshel); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    call Nshel_num_basis(nshel, num_basis); check_err()
    allocate(ao(num_basis, 1))

    write(ifile_out,'("i,j,val")') 
    read(ifile_in,*)
    ir = 1
    do
       read(ifile_in, *, end=100) rs(1,1), rs(1,2), rs(1,3)
       call Nshel_ao_at(nshel, rs, n_dr, ao); check_err()
       do mu = 1, num_basis
          write(ifile_out,'(i0,",",i0,",",f20.10)') mu, ir, ao(mu,1)
       end do
       ir = ir+1
    end do
    
100 continue

    close(ifile_in)
    close(ifile_out)        
    
  end subroutine AOGrid_run
  subroutine AOAOGrid_run
    use Mod_fjson
    use Mod_ArgParser
    type(Obj_Nshel) :: nshel
    character(100) :: fn_nshel, fn_in, fn_out, str_type
    double precision :: x, y
    double precision, allocatable :: zs(:)
    integer :: ir, num_basis, mu, nu
    integer :: nj_dz, nk_dz, nr
    double precision, allocatable :: aoao(:,:,:)
    integer, parameter :: ifile_in = 192
    integer, parameter :: ifile_out = 193
    
    call arg_parse_s("-t", str_type); check_err()
    call arg_parse_s("-i", fn_in); check_err()
    call arg_parse_s("-o", fn_out); check_err()
    call arg_parse_s("-nshel", fn_nshel); check_err()

    select case(str_type)
    case("ii00")
       nj_dz = 0
       nk_dz = 0
    case("ii01")
       nj_dz = 0
       nk_dz = 1
    case default
       throw_err("unsupported str_type", 1)
    end select
    write(*,*) "x: integrate"
    write(*,*) "y: integrate"
    write(*,*) "z:", nj_dz, nk_dz

    call open_r(ifile_in, fn_in); check_err()
    call open_w(ifile_out, fn_out); check_err()
    
    write(*,*) "nshel json file:", trim(fn_nshel)
    call Nshel_new_file(nshel, fn_nshel); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    call Nshel_num_basis(nshel, num_basis); check_err()    

    nr = 0
    read(ifile_in, *)
    do
       read(ifile_in, *, end=101)
       nr = nr +1
    end do
101 continue

    allocate(zs(nr))
    allocate(aoao(num_basis, num_basis, nr))
    
    rewind(ifile_in)
    read(ifile_in, *)
    do ir = 1, nr
       read(ifile_in, *) x, y, zs(ir)
    end do

    call Nshel_aoao_at_intxy(nshel, zs(:), nj_dz, nk_dz, aoao(:,:,:)); check_err()

    write(ifile_out,'("i,j,k,val")') 
    do ir = 1, nr
       do mu = 1, num_basis
          do nu = 1, num_basis
             write(ifile_out,'(i0,",",i0,",",i0,",",f20.10)') mu, nu,ir, aoao(mu,nu,1)
          end do
       end do
    end do
    
    close(ifile_in)
    close(ifile_out)        
    
  end subroutine AOAOGrid_run
  ! -- AO matrix --
  subroutine AOMat_run
    use Mod_fjson
    use Mod_ArgParser
    use Mod_math, only : dump_dmat
    type(Obj_Nshel) :: nshel
    character(100) :: fn_nshel, out_dir, str_type
    integer :: num_basis
    double precision, allocatable :: m(:,:), m3(:,:,:)
    integer :: ifile = 192
    
    call arg_parse_s("-t", str_type); check_err()
    call arg_parse_s("-out_dir", out_dir); check_err()
    call arg_parse_s("-nshel", fn_nshel); check_err()    
    
    write(*,*) "nshel json file:", trim(fn_nshel)
    call Nshel_new_file(nshel, fn_nshel); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    call Nshel_num_basis(nshel, num_basis); check_err()

    allocate(m(num_basis, num_basis), m3(3, num_basis, num_basis))

    select case(str_type)
    case("dw2")
       call Nshel_dw2(nshel, m3); check_err()

       call open_w(ifile, trim(out_dir)//"/ao_dx2.csv"); check_err()       
       call dump_dmat(m3(1,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call open_w(ifile, trim(out_dir)//"/ao_dy2.csv"); check_err()       
       call dump_dmat(m3(2,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call open_w(ifile, trim(out_dir)//"/ao_dz2.csv"); check_err()       
       call dump_dmat(m3(3,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
    case("dw")
       call Nshel_dw(nshel, m3); check_err()

       call open_w(ifile, trim(out_dir)//"/ao_dx.csv"); check_err()       
       call dump_dmat(m3(1,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call open_w(ifile, trim(out_dir)//"/ao_dy.csv"); check_err()       
       call dump_dmat(m3(2,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call open_w(ifile, trim(out_dir)//"/ao_dz.csv"); check_err()       
       call dump_dmat(m3(3,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

    case("r")
       call Nshel_r(nshel, m3); check_err()
       
       call open_w(ifile, trim(out_dir)//"/ao_rx.csv"); check_err()       
       call dump_dmat(m3(1,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_ry.csv"); check_err()
       call dump_dmat(m3(2,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_rz.csv"); check_err()       
       call dump_dmat(m3(3,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
    case("s")
       call Nshel_s(nshel, m); check_err()

       call open_w(ifile, trim(out_dir)//"/ao_s.csv"); check_err()       
       call dump_dmat(m(:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1

    case("stvrd")
       call Nshel_s(nshel, m); check_err()
       call open_w(ifile, trim(out_dir)//"/ao_s.csv"); check_err()       
       call dump_dmat(m, ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call Nshel_t(nshel, m); check_err()
       call open_w(ifile, trim(out_dir)//"/ao_t.csv"); check_err()       
       call dump_dmat(m, ifile); check_err()
       close(ifile)
       ifile = ifile + 1

       call Nshel_v(nshel, m); check_err()
       call open_w(ifile, trim(out_dir)//"/ao_v.csv"); check_err()       
       call dump_dmat(m, ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call Nshel_r(nshel, m3); check_err()
       call open_w(ifile, trim(out_dir)//"/ao_rx.csv"); check_err()       
       call dump_dmat(m3(1,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_ry.csv"); check_err()
       call dump_dmat(m3(2,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_rz.csv"); check_err()       
       call dump_dmat(m3(3,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call Nshel_dw(nshel, m3); check_err()
       call open_w(ifile, trim(out_dir)//"/ao_dx.csv"); check_err()
       call dump_dmat(m3(1,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_dy.csv"); check_err()
       call dump_dmat(m3(2,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
       
       call open_w(ifile, trim(out_dir)//"/ao_dz.csv"); check_err()
       call dump_dmat(m3(3,:,:), ifile); check_err()
       close(ifile)
       ifile = ifile + 1
    case default
       throw_err("unsupported", 1)
    end select
        
  end subroutine AOMat_run
  subroutine print_help()
    write(0,*) "argument is necessary"
    write(0,*) "endvr runtype nshel.json rs.dat mat.dat"
  end subroutine print_help
end module Mod_Enint

program main
  use Mod_Enint
  call Enint_run
  if(get_err().ne.0) then
     write(0,*) "Error on calculation. Stop enint program..."
     call exit(1)
  end if
end program main
