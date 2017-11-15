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
  ! -- produce Grid --
  subroutine Grid_run
    use Mod_ArgParser
    use Mod_StrUtil
    character(100) :: line, fn_out
    integer ix, iy, iz, nx, ny, nz
    double precision, allocatable :: xs(:), ys(:), zs(:)
    integer, parameter :: ifile = 411
    
    call arg_parse_s("-x", line); check_err()
    call str2vec(line, nx, xs, .false.); check_err()
    allocate(xs(nx))
    call str2vec(line, nx, xs, .true.); check_err()

    call arg_parse_s("-y", line); check_err()
    call str2vec(line, ny, ys, .false.); check_err()
    allocate(ys(ny))
    call str2vec(line, ny, ys, .true.); check_err()

    call arg_parse_s("-z", line); check_err()
    call str2vec(line, nz, zs, .false.); check_err()
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
  ! -- AO grid representation --
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
    end select        

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
  ! -- AO matrix --
  subroutine AOMat_run
    use Mod_fjson
    use Mod_ArgParser
    type(Obj_Nshel) :: nshel
    character(100) :: fn_nshel, fn_in, out_dir, str_type
    double precision :: rs(1,3)
    integer :: ir, num_basis, mu
    integer :: n_dr(3)
    double precision, allocatable :: m(:,:)
    integer :: ifile = 192
    
    call arg_parse_s("-t", str_type); check_err()
    call arg_parse_s("-out_dir", out_dir); check_err()
    call arg_parse_s("-nshel", fn_nshel); check_err()    
    
    write(*,*) "nshel json file:", trim(fn_nshel)
    call Nshel_new_file(nshel, fn_nshel); check_err()
    call Nshel_setup(nshel, .true.); check_err()
    call Nshel_num_basis(nshel, num_basis); check_err()

    allocate(m(num_basis, num_basis))

    select case(str_type)
    case("stvrd")
       call Nshel_s(nshel, m); check_err()
       call open_w(ifile, out_dir//"/ao_s.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_t(nshel, m); check_err()
       call open_w(ifile, out_dir//"/ao_t.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_v(nshel, m); check_err()
       call open_w(ifile, out_dir//"/ao_v.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_r(nshel, 1, m); check_err()
       call open_w(ifile, out_dir//"/ao_rx.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_r(nshel, 2, m); check_err()
       call open_w(ifile, out_dir//"/ao_ry.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_r(nshel, 3, m); check_err()
       call open_w(ifile, out_dir//"/ao_rz.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()                     

       call Nshel_dw(nshel, 1, m); check_err()
       call open_w(ifile, out_dir//"/ao_dx.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_dw(nshel, 2, m); check_err()
       call open_w(ifile, out_dir//"/ao_dy.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()

       call Nshel_dw(nshel, 3, m); check_err()
       call open_w(ifile, out_dir//"/ao_dz.csv"); check_err()
       ifile = ifile + 1
       call dump_dmat(m, ifile); check_err()       

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
