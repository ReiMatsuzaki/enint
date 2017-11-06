#include "macros.fpp"
#include "macros_utest.fpp"
module Mod_TestArgParser
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_ArgParser
  implicit none
contains
  subroutine TestArgParser_run

 !   call Utest_sub_begin("first")
 !   call test_first
 !   call Utest_sub_end

    call Utest_sub_begin("second")
    call test_second
    call Utest_sub_end    
    
  end subroutine TestArgParser_run
  subroutine test_first()
    type(object) :: o
    call ArgParser_new
    call ArgParser_add("s", TYPE_S, 1)
    call ArgParser_add("i", TYPE_I, 1)
    call ArgParser_add("x", TYPE_D, 1)
    call ArgParser_add("y", TYPE_D, 1)
    call ArgParser_parse(o)

    call object_dump(o)
    
    call ArgParser_delete
  end subroutine test_first
  subroutine test_second()
    integer :: i
    double precision :: x, y(3)
    character(100) :: s
    
    call arg_parse_i("-i", i); check_err()
    call expect_eq(523, i); check_err()

    call arg_parse_s("-s", s); check_err()
    call expect_eq("s1s2", s); check_err()

    call arg_parse_d("-x", x); check_err()
    call expect_eq(10.1d0, x); check_err()

    call arg_parse_dvec("-y", y); check_err()
    call expect_eq(10.2d0, y(1)); check_err()
    call expect_eq(10.3d0, y(2)); check_err()
    call expect_eq(10.4d0, y(3)); check_err()
    
  end subroutine test_second
end module Mod_TestArgParser

program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestArgParser

  call utest_new
  call ErrHandle_new
  
  call TestArgParser_run

  call ErrHandle_delete
  call Timer_result
  call utest_delete
  
end program main
