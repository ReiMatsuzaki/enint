#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_TestDet
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_DetString
  implicit none
contains
  subroutine TestDet_run
    
    call Utest_sub_begin("build")
    call test_build()
    call Utest_sub_end()
    
  end subroutine TestDet_run
  subroutine test_build

    !n\N|  0 1 2
    !---+------------
    ! 0 |  1
    !   |  |\
    ! 1 |  2 3
    !   |  |\|\
    ! 2 |  4 5 6
    !   |   \|\|
    ! 3 |    7 8
    !   |     \|
    ! 4 |      9
    type(Obj_DetString) :: str
    call DetString_new(str, 2, 4)
    call DetString_dump(str)
    
    call expect_true(str%o(0,0))
    call expect_false(str%o(0,1))
    call expect_false(str%o(0,2))

    call expect_true(str%o(1,0))
    call expect_true(str%o(1,1))
    call expect_false(str%o(1,2))

    call expect_true(str%o(2,0))
    call expect_true(str%o(2,1))
    call expect_true(str%o(2,2))

    call expect_false(str%o(3,0))
    call expect_true(str%o(3,1))
    call expect_true(str%o(3,2))

    call expect_false(str%o(4,0))
    call expect_false(str%o(4,1))
    call expect_true(str%o(4,2))

    call expect_eq(3, str%w(3, 1))
    call expect_eq(3, str%y(4,2))

  end subroutine test_build
end module Mod_TestDet

program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestDet

  call utest_new
  call ErrHandle_new
  
  call TestDet_run

  call ErrHandle_delete
  call Timer_result
  call utest_delete
  
end program main
