#include "macros.fpp"
#include "macros_utest.fpp"
module Mod_TestStrUtil
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_StrUtil
  implicit none
contains
  subroutine TestStrUtil_run
    call Utest_sub_begin("str_split")
    call test_str_split
    call Utest_sub_end

    call Utest_sub_begin("str2vec")
    call test_str_split
    call Utest_sub_end    
  end subroutine TestStrUtil_run
  subroutine test_str_split
    character(20) :: line
    character(20) :: lines(10)
    integer n
    line = "x1,x2,3.14"
    call str_split(line, ",", n, lines(:))
    call expect_eq(3, n)
    call expect_eq("x1", lines(1))
    call expect_eq("x2", lines(2))
    call expect_eq("3.14", lines(3))

    line = "momo"
    call str_split(line, "n", n,lines(:))
    call expect_eq(1, n)
    call expect_eq("momo", lines(1))
    
  end subroutine test_str_split
  subroutine test_str2vec
    character(100) :: line
    double precision :: xs(5)
    integer n
    line = "linspace:1.0:3.0:3"
    call str2vec(line, n, xs); check_err()
    call expect_eq(n,3)
    call expect_eq(1.0d0, xs(1))
    call expect_eq(2.0d0, xs(2))
    call expect_eq(3.0d0, xs(3))

    line = "scalar:1.12"
    call str2vec(line, n, xs); check_err()
    call expect_eq(n,1)
    call expect_eq(1.12d0, xs(1))
    
  end subroutine test_str2vec
end module Mod_TestStrUtil

program main
  use Mod_ErrHandle
  use Mod_Utest
  use Mod_UtestCheck
  use Mod_TestStrUtil

  call utest_new
  call ErrHandle_new
  
  call TestStrUtil_run

  call ErrHandle_delete
  call Timer_result
  call utest_delete
  
end program main
