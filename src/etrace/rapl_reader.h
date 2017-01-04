/*
 * A simple wrapper code for rapl_msr
 * 
 * Kaz Yoshii <ky@anl.gov>
 *
 * NOTE: this is not thread safe.
 */

#ifndef __DEFINE_READ_RAPL_H_DEFINED__
#define __DEFINE_READ_RAPL_H_DEFINED__

#define RAPL_READER_VERSION  "0.2"

#define USE_MSR_SAFE

#include <stdint.h>

#include "msr-index.h"

/* 
   See  Intel 64 and IA-32 Architectures
   Software Developer’s Manual
   Volume 3B: System Programming Guide, Part 2
   p21-

   NOTE:
   MSR_RAPL_POWER_UNIT provides common information for all domains
   MSR_(PKG|DRAM)_ENERGY_STATUS is updated every ~1ms
*/

typedef  struct {
  uint32_t  pkg, pp0, pp1, mem; /* 32-bit wide */
} rapl_domain_energy_t;

typedef  struct {
  double  pkg, pp0, pp1, mem;
} rapl_domain_energy_total_t;


typedef  struct {
  int fd, coreid;

  uint64_t  pkg_power_info;
  uint64_t  dram_power_info;

  uint64_t  pkg_power_limit;
  uint64_t  pkg_power_limit_new;
  uint64_t  dram_power_limit;
  uint64_t  pp0_power_limit;
  uint64_t  pp1_power_limit;

  uint64_t  pp0_policy;
  uint64_t  pp1_policy;

  rapl_domain_energy_t  e[2];  /* enery snapshot two element circular buffer */
  double timestamp[2]; /* timestamp associated with the snapshot */
  int pos; /* position of e[].  0 or 1 */

  int count; /* snapshot count */

  rapl_domain_energy_t  e_start; /* enery snapshot at init */

  rapl_domain_energy_total_t  e_total; /* summed up at internal sampling */

  double timestamp_start; 
} rapl_pkg_t;

#define MAX_RAPL_PKG (16)   /* # of sockets */

typedef struct {
  int model;

  double power_units,  energy_units, time_units;

  int read_dram, read_pp1;

  rapl_pkg_t pkg[MAX_RAPL_PKG];

  int npkg;
} rapl_reader_t;




/* rapl_reader_init() return 0 on success, otherwise
   it returns non-zero */
extern int   rapl_reader_init(void);
extern void  rapl_reader_snap(void);
extern void  rapl_reader_start(void);
extern void  rapl_reader_stop(void);
extern void  rapl_reader_report_abs(void);
extern void  rapl_reader_report_delta(void);
extern void  rapl_reader_finalize(void);

extern int   rapl_reader_domain_pp1(void);
extern int   rapl_reader_domain_dram(void);

extern void rapl_reader_print_info_header(FILE* fp);
extern void rapl_reader_print_info_footer(FILE* fp);

extern void rapl_set_power_limit(double plim1, double plim2);


/*
  pkg, pp0, pp1 and mem are an array of double. e.g. pp0[1] means
  power consumption of pp0 in the second package(or socket).  power
  consumption in Watt (J/s)
 */
extern int rapl_reader_get_power(double *timestamp, double *pkg, double *pp0, 
				 double *pp1, double *mem );

extern const rapl_reader_t*  rapl_reader_get_struct(void);

extern int rapl_reader_get_energy(double *_delta, double *pkg, double *pp0, 
				  double *pp1, double *mem );

#endif
