#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#ifdef __cplusplus
extern "C" {
#endif

#define __IN__
#define __OUT__ 
#define __IN_OUT__


#define _D2_DOUBLE
#define _VERBOSE_OUTPUT

#ifdef _D2_DOUBLE
#define SCALAR double
#define SCALAR_STDIO_TYPE ("%lf ")
#elif defined _D2_SINGLE
#define SCALAR float
#define SCALAR_STDIO_TYPE ("%f ")
#endif

  extern int world_rank; // rank of processor

#ifdef _VERBOSE_OUTPUT
#include <stdarg.h>
#include <stdio.h>
  static inline void VPRINTF(const char *format, ...) {
    va_list args;
#ifdef __USE_MPI__
    if (world_rank == 0) {
#endif
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
#ifdef __USE_MPI__
    }
#endif
  }

  static inline void VFLUSH() {fflush(stdout);}
#else
  static inline void VPRINTF(const char *format, ...) {}
  static inline void VFLUSH() {}
#endif


#ifdef  _D2_DOUBLE
#define _D2_SCALAR          double
#define _D2_FUNC(x)         _d ## x
#define _D2_CBLAS_FUNC(x)   cblas_d ## x
  //#define _D2_LAPACKE_FUNC(x) d ## x
#elif defined  _D2_SINGLE
#define _D2_SCALAR          float
#define _D2_FUNC(x)         _s ## x
#define _D2_CBLAS_FUNC(x)   cblas_s ## x
  //#define _D2_LAPACKE_FUNC(x) s ## x
#endif


// Timing, count in nano seconds.
#include <time.h>

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
#include <mach/clock.h>
#include <mach/mach.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
inline int clock_gettime(int clk_id, struct timespec* ts) {
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
  return 0;
}
#endif

#define BILLION  1000000000L
static struct timespec nstart, nend;
static inline void nclock_start() {clock_gettime(CLOCK_MONOTONIC, &nstart);}
static inline double nclock_end() {clock_gettime(CLOCK_MONOTONIC, &nend);     
  return (double) ( nend.tv_sec - nstart.tv_sec ) + (double) ( nend.tv_nsec - nstart.tv_nsec ) / (double) BILLION ;
}

#ifdef __cplusplus
}
#endif

#endif /* _GLOBAL_H_ */
