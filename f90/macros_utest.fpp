#ifndef MACROS_UTEST_FPP

      
!     Unit test
#if defined GOOD_MACRO
#define expect_near(a,b,eps) utest_check_near_0((a),(b),(eps),(__FILE__),(__LINE__),.FALSE.)
#define assert_near(a,b,eps) utest_check_near_0((a),(b),(eps),(__FILE__),(__LINE__),.TRUE.)

#define expect_eq(a,b) utest_check_eq_0((a),(b),(__FILE__),(__LINE__),.FALSE.)
#define assert_eq(a,b) utest_check_eq_0((a),(b),(__FILE__),(__LINE__),.TRUE.)

#define expect_true(a)  utest_check_true((a), (#a), (__FILE__), (__LINE__), .FALSE.)
#define assert_true(a)  utest_check_true((a), (#a), (__FILE__), (__LINE__), .TRUE.)
#define expect_false(a)  utest_check_false((a), (#a), (__FILE__), (__LINE__), .FALSE.)
#define assert_false(a)  utest_check_false((a), (#a), (__FILE__), (__LINE__), .TRUE.)
      
#else
#define expect_near(a,b,eps) utest_check_near_0((a),(b),(eps),(__FILE__),(__LINE__),.FALSE.)
#define assert_near(a,b,eps) utest_check_near_0((a),(b),(eps),(__FILE__),(__LINE__),.TRUE.)

#define expect_eq(a,b) utest_check_eq_0((a),(b),(__FILE__),(__LINE__),.FALSE.)
#define assert_eq(a,b) utest_check_eq_0((a),(b),(__FILE__),(__LINE__),.TRUE.)

#define expect_true(a)  utest_check_true_0((a), (__FILE__), (__LINE__), .FALSE.)
#define assert_true(a)  utest_check_true_0((a), (__FILE__), (__LINE__), .TRUE.)
#define expect_false(a)  utest_check_false_0((a), (__FILE__), (__LINE__), .FALSE.)
#define assert_false(a)  utest_check_false_0((a), (__FILE__), (__LINE__), .TRUE.)

#endif

#endif
