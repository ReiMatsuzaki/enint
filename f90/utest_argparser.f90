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

    call Utest_sub_begin("first")
    call test_first
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
