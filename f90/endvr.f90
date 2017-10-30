#include "macros.fpp"

module Mod_Endvr
  ! compute electron nuclear interaction matrix when
  ! electron part is represented by GTO and
  ! nuclear part is represented by DVR
  use Mod_ErrHandle
  use Mod_Nshel
  use Mod_Fjson
  use Mod_Timer
  implicit none
  type(Obj_Nshel) :: nshel_
  double precision :: q_ ! charge of nuclear
  double precision, allocatable :: rs_(:,:) ! rs_(n,:) : n th path
  
contains
  subroutine Endvr_run
    integer, parameter :: iin=123, iout=124
    character(100) :: fn
    double precision, allocatable :: v(:,:)
    integer :: num, a, i, j, iargc

    call ErrHandle_new()
    call Timer_new("Endvr", .false.); check_err()

    write(*,*) " ==== Endvr ===="

    if(iargc()<3) then
       begin_err(1)
       call print_help
       end_err()
    end if

    call getarg(1, fn)
    write(*,*) "nshel json file:", fn
    call Nshel_new_file(nshel_, fn); check_err()
    call Nucs_delete(nshel_%nucs); check_err()
    call nucs_new(nshel_%nucs, 1); check_err()    
    nshel_%nucs%zs(1) = 1.0d0
    call Nshel_setup(nshel_); check_err()
    num = nshel_%nbasis
    allocate(v(num, num))

    call getarg(2, fn)
    write(*,*) "input rs file:", fn
    call open_r(iin, fn); check_err()

    call getarg(3, fn)
    write(*,*) "output file: ", fn
    call open_w(iout, fn); check_err()

    ! - start read -
    read(iin, *)
    write(iout,*) "a,i,j,val"
    a = 0
    do
       a = a+1
       
       read(iin, *, end=100) nshel_%nucs%ws(1,:)
       
       call Nshel_v(nshel_, v(:,:)); check_err()

       do i = 1, num
          do j = 1, num
             write(iout,'(i0,",",i0,",",i0,f20.10)') a, i, j, v(i,j)
          end do
       end do
    end do

100 continue
    
    close(iin)
    close(iout)
    call Nshel_delete(nshel_); check_err()
    call Timer_result()
    call Timer_delete()
    call ErrHandle_delete()
    
  end subroutine Endvr_run
  subroutine print_help()
    write(0,*) "argument is necessary"
    write(0,*) "endvr nshel.json rs.dat mat.dat"
  end subroutine print_help
end module Mod_Endvr

program main
  use Mod_Endvr
  call Endvr_run
end program main
