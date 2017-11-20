#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_TestDRT
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_DRT
  implicit none
contains
  subroutine TestDRT_run
    
    call Utest_sub_begin("build")
    call test_build()
    call Utest_sub_end()

    call Utest_sub_begin("wk")
    call test_wk()
    call Utest_sub_end()    
    
  end subroutine TestDRT_run
  subroutine test_build
    type(Obj_DRT) :: drt
    call DRT_new(drt, 4,0,4)

    call expect_eq(20, drt%xj(1))
    call expect_eq(6,  drt%xj(2))
    call expect_eq(8,  drt%xj(3))
    call expect_eq(6,  drt%xj(4))
    call expect_eq(2,  drt%xj(9))
    
    !  call DRT_dump(drt, .true.)
    !  call DRT_dump(drt, .false.)

    call expect_eq(4, ilevel(drt, 1))
    call expect_eq(3, ilevel(drt, 2))
    call expect_eq(2, ilevel(drt, 7))
    call expect_eq(1, ilevel(drt, 13))
    call expect_eq(0, ilevel(drt, 14))
    
  end subroutine test_build
  subroutine test_wk
    integer, parameter :: norbs = 4
    type(Obj_DRT) :: drt
    integer :: ds(norbs)
    integer lwk, js(0:norbs)
    logical :: findq

    call DRT_new(drt, 4,0,4)
    
    ds(:) = NANVAL
    ds(3) = 0
    ds(2) = 3
    call expect_true(correct_wk(drt, 3,12, ds))

    ds(4) = 2
    ds(3) = 1
    ds(2) = 2
    ds(1) = 1
    call nodes_wk(drt, ds, js)
    call expect_eq(js(4), 1)
    call expect_eq(js(3), 3)
    call expect_eq(js(2), 7)
    call expect_eq(js(1), 12)
    call expect_eq(js(0), 14)
    
    ds(4) = 0
    ds(3) = 0
    ds(2) = 3
    ds(1) = 3
    do lwk = 1, num_wk(drt)
       call expect_true(correct_wk(drt, 1,drt%nj, ds))
       call expect_eq(lwk, index_wk(drt, ds))
       call next_wk(drt, ds, findq)
       if(lwk.ne.num_wk(drt)) then
          call expect_true(findq)
       else
          call expect_false(findq)
       end if          
    end do

  end subroutine test_wk
end module Mod_TestDRT

program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestDRT

  call utest_new
  call ErrHandle_new
  
  call TestDRT_run

  call ErrHandle_delete
  call Timer_result
  call utest_delete
  
end program main
