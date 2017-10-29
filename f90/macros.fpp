#ifndef MACROS_FPP
      
!     Error handling
#define throw_err(msg, ierr0) call ErrHandle_err(msg, ierr0, (__FILE__), (__LINE__)); return
#define throw_err_stop(msg, ierr0) call ErrHandle_err(msg, ierr0, (__FILE__), (__LINE__)); stop
#define begin_err(ierr0) call ErrHandle_err("", ierr0, (__FILE__), (__LINE__));
#define end_err() return
#define check_err() if(ierr .ne. 0) begin_err(ierr); if(ierr .ne. 0) return
#define check_err_stop() if(ierr.ne.0) begin_err(ierr); if(ierr.ne.0) stop
#define get_err() ierr
  
#endif
