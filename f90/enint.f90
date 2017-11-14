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
    write(*,*) "runtype:", runtype

    select case(runtype)
    case("na_rs")
       call na_rs_run; check_err()
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
  ! -- AO grid representation --
  subroutine AOGrid_run
    use Mod_fjson
    character(100) :: fn_nshel, line

    if(iargc()<4) then
       begin_err(1)
       call print_help
       end_err()
    end if

    call getarg(2, fn_nshel)
    write(*,*) "nshel json file:", fn_nshel
    call getarg(3, line)
    
  end subroutine AOGrid_run
  subroutine print_help()
    write(0,*) "argument is necessary"
    write(0,*) "endvr runtype nshel.json rs.dat mat.dat"
  end subroutine print_help
end module Mod_Enint

program main
  use Mod_Enint
  call Enint_run
end program main
