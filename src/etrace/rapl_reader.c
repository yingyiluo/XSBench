/*
  A simple wrapper code for rapl_msr

  requires the LLNL msr_safe driver
  
  Kazutomo Yoshii <ky@anl.gov>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>

#include "rapl_reader.h"


double gettimesec(void)
{
	struct timeval tv;

	gettimeofday(&tv, 0);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000.0/1000.0;
}


/* from parse_cpuinfo_pkg_model.c */
extern int parse_cpuinfo_pkg_model(int *cores, int ncores,int *model);

static rapl_reader_t rapl;

#if defined(USE_MSR_SAFE)
static int  _open_msr(int coreid)
{
  char fn[1024];
  int fd;

  snprintf(fn,sizeof(fn),"/dev/cpu/%d/msr_safe",coreid);
  fd = open(fn, O_RDWR);
  if( fd<0 ) {
    printf("\n");
    printf("Failed to open %s. Is the msr safe driver installed? If so, check the permission.\n", fn);
    printf("\n");
    return -1;
  }
  return fd;
}
#else
static int  _open_msr(int coreid)
{
  char fn[1024];
  int fd;

  snprintf(fn,sizeof(fn),"/dev/cpu/%d/msr",coreid);
  fd = open(fn, O_RDWR);
  if( fd<0 ) {
    printf("\n");
    printf("failed to open %s!\n", fn);
    printf("\n");
    printf("[Suggestions]\n");
    printf("\n");
    printf("sudo chmod o+rw %s\n", fn);
    printf("\n");
    printf("You might also need to run setcap.\n");
    printf("\n");
    printf("If you use the etrace command:\n");
    printf("\n");
    printf("sudo setcap cap_sys_rawio=ep etrace\n");
    printf("\n");
    printf("If you call rapl_reader library functions from your application:\n");
    printf("\n");
    printf("sudo setcap cap_sys_rawio=ep your_binary\n");
    printf("\n");
    printf("\n");
    return -1;
  }
  return fd;
}
#endif


static uint64_t _read_msr(int fd, int offset) 
{
  uint64_t data;

  if( pread(fd, &data, 8, offset) != 8 ) {
	  printf("pread at %x\n", offset);
	  return -1;
  }
  return data;
}

static uint64_t _write_msr(int fd, int offset, uint64_t data) 
{
  if( pwrite(fd, &data, 8, offset) != 8 ) {
    perror("pwrite");
    return -1;
  }
  return data;
}



int   rapl_reader_domain_pp1(void)
{
  return rapl.read_pp1;
}

int   rapl_reader_domain_dram(void)
{
  return rapl.read_dram;
}

/*
  Extract (h-l+1) bits from val and return l-bit shifted value
*/

inline uint64_t extractbits(uint64_t val, int l, int h )
{
  int w= h - l+ 1;
  uint64_t m;
  if(w<1) return 0;

  val = val >> l;
  m = (1UL<<w)-1;
  val = val & m;

  return val;
}

/*
  Insert newval into (h-l+1) bits of val and return updated val
*/
inline uint64_t insertbits(uint64_t val, int l, int h,uint64_t newval )
{
  int w= h - l+ 1;
  uint64_t m;
  if(w<1) return 0;

  m = (1UL<<w)-1;
  newval = newval & m;

  m = (1UL<<(h+1))-1 - ((1UL<<l)-1);
  m = ~m;

  val = (val&m) | (newval<<l);

  return val;
}
/* ad hoc */
static uint64_t prev_aperf[16];
static uint64_t prev_mperf[16];

void  print_adhoc(FILE* fp)
{
  int i;
  uint64_t t1,t2;

  for(i=0;i<rapl.npkg;i++) {
    t1 = _read_msr(rapl.pkg[i].fd, MSR_IA32_MPERF);
    t2 = _read_msr(rapl.pkg[i].fd, MSR_IA32_APERF);
    if( prev_aperf[i] == 0 ) {
      fprintf(fp,"0.0 ");
    } else {
      fprintf(fp,"%lf ", 
	      (double)(t1-prev_aperf[i])/
	      (double)(t2-prev_mperf[i]) );
    }
    prev_aperf[i] = t1;
    prev_mperf[i] = t2;
  }


  /*
  fprintf(fp,"PKG:");
  for(i=0;i<rapl.npkg;i++) {
    tmpval = _read_msr(rapl.pkg[i].fd, MSR_PKG_PERF_STATUS) & 0xffffffff;
    fprintf(fp,"%d ",(int)tmpval);
  }
  fprintf(fp,"PP0:");
  for(i=0;i<rapl.npkg;i++) {
    uint64_t tmpval;
    tmpval = _read_msr(rapl.pkg[i].fd, MSR_PP0_PERF_STATUS) & 0xffffffff;
    fprintf(fp,"%d ",(int)tmpval);
  }
  fprintf(fp,"DRAM:");
  for(i=0;i<rapl.npkg;i++) {
    uint64_t tmpval;
    tmpval = _read_msr(rapl.pkg[i].fd, MSR_DRAM_PERF_STATUS) & 0xffffffff;
    fprintf(fp,"%d ",(int)tmpval);
  }
  fprintf(fp,"PT:");
  for(i=0;i<rapl.npkg;i++) {
    uint64_t tmpval;
    tmpval = _read_msr(rapl.pkg[i].fd, MSR_IA32_PACKAGE_THERM_STATUS);
    fprintf(fp,"%lx ",tmpval);
  }
  */
}


static double _get_power_info_thermal_spec(uint64_t pi)
{
  return rapl.power_units*(double)extractbits(pi, 0,14);
}
static double _get_power_info_min_power(uint64_t pi)
{
  return rapl.power_units*(double)extractbits(pi,16,30);
}
static double _get_power_info_max_power(uint64_t pi)
{
  return rapl.power_units*(double)extractbits(pi,32,46);
}
static double _get_power_info_time_window(uint64_t pi)
{
  return rapl.time_units*(double)extractbits(pi,48,53);
}

static double _get_power_limit_power1(uint64_t pi)
{
  return rapl.power_units*(double)extractbits(pi, 0,14);
}
static uint64_t   _set_power_limit_power1(uint64_t pi, double newlimit)
{
  uint64_t val = (uint64_t)newlimit/rapl.power_units;
  return insertbits(pi, 0,14,val);
}
static int _get_power_limit_enable1(uint64_t pi)
{
  return extractbits(pi,15,15);
}
static int _get_power_limit_clamp1(uint64_t pi)
{
  return extractbits(pi,16,16);
}
static uint64_t _set_power_limit_clamp1(uint64_t pi,int newval)
{
  return insertbits(pi,16,16,newval);
}
static double _get_power_limit_time1(uint64_t pi)
{
  return rapl.time_units*(double)extractbits(pi,17,23);
}



static double _get_power_limit_power2(uint64_t pi)
{
  return rapl.power_units*(double)extractbits(pi,32,46);
}
static uint64_t   _set_power_limit_power2(uint64_t pi, double newlimit)
{
  uint64_t val = (uint64_t)newlimit/rapl.power_units;
  return insertbits(pi,32,46,val);
}

static int _get_power_limit_enable2(uint64_t pi)
{
  return extractbits(pi,47,47);
}
static int _get_power_limit_clamp2(uint64_t pi)
{
  return extractbits(pi,48,48);
}
static uint64_t _set_power_limit_clamp2(uint64_t pi,int newval)
{
  return insertbits(pi,48,48,newval);
}
static double _get_power_limit_time2(uint64_t pi)
{
  return rapl.time_units*(double)extractbits(pi,49,55);
}
static int _get_power_limit_lock(uint64_t pi)
{
  return extractbits(pi,63,63);
}

void rapl_set_power_limit(double plim1, double plim2)
{
  int i;

  printf("plim1=%lf plim2=%lf\n",plim1,plim2);
  for(i=0;i<rapl.npkg;i++) {
    uint64_t tmpval = rapl.pkg[i].pkg_power_limit;

    tmpval = _set_power_limit_power1(tmpval, plim1);
    tmpval = _set_power_limit_clamp1(tmpval, 1 );
    tmpval = _set_power_limit_power2(tmpval, plim2);
    tmpval = _set_power_limit_clamp2(tmpval, 1 );

    _write_msr(rapl.pkg[i].fd, MSR_PKG_POWER_LIMIT, tmpval );
    rapl.pkg[i].pkg_power_limit_new = tmpval;

    /* debug */
    if(0) {
      tmpval = _read_msr(rapl.pkg[i].fd, MSR_PKG_POWER_LIMIT);

      fprintf(stdout,"@@ SOCKET%d_POWER_LIMIT: LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LOCK=%d\n",i,
	    _get_power_limit_power1(tmpval),
	    _get_power_limit_enable1(tmpval),
	    _get_power_limit_clamp1(tmpval),
	    _get_power_limit_time1(tmpval),
	    _get_power_limit_power2(tmpval),
	    _get_power_limit_enable2(tmpval),
	    _get_power_limit_clamp2(tmpval),
	    _get_power_limit_time2(tmpval),
	    _get_power_limit_lock(tmpval)
	    );
    }
  }
}


void rapl_reset_power_limit(void)
{
  int i;
  for(i=0;i<rapl.npkg;i++) {
    uint64_t tmpval = rapl.pkg[i].pkg_power_limit;
    _write_msr(rapl.pkg[i].fd, MSR_PKG_POWER_LIMIT, tmpval );
  }
}



int rapl_reader_init(void)
{
  int cores[MAX_RAPL_PKG]; /* uniq core per pkg */
  int i;
  uint64_t val;
  int fd;

  memset(&rapl, 0, sizeof(rapl));

  //#define DEBUGTMP
#ifdef DEBUGTMP
  {
#define NCPUS (32)
    int fds[NCPUS];
    int j;
    uint64_t  pkgval, pp0val, pp1val, dramval;
    for(i=0; i<NCPUS; i++ ) {
      fds[i] = _open_msr(i);
    }
    printf("%lx\n", _read_msr(fds[0] ,MSR_RAPL_POWER_UNIT));

    for(j=0; j<10; j++ ) {
      for(i=0; i<NCPUS; i++ ) {
	pkgval=_read_msr( fds[i], MSR_PKG_ENERGY_STATUS);
	pp0val=_read_msr( fds[i], MSR_PP0_ENERGY_STATUS);
	//pp1val=_read_msr( fds[i], MSR_PP1_ENERGY_STATUS);
	dramval=_read_msr(fds[i], MSR_DRAM_ENERGY_STATUS);
	printf("%2d %lu %lu %lu\n", i, pkgval, pp0val,  dramval);
      }
      sleep(1);
      puts("");
    }

    for(i=0; i<32; i++ ) {
      close(fds[i]); 
    }

  }
#endif


  rapl.npkg = parse_cpuinfo_pkg_model(cores, MAX_RAPL_PKG, &rapl.model );

  for(i=0;i<rapl.npkg;i++) {
    rapl.pkg[i].coreid = cores[i]; /* unique coreid per package */
    fd = _open_msr( rapl.pkg[i].coreid );
    if( fd < 0 ) {
	    return -1;
    }
    rapl.pkg[i].fd = fd;
  }

  if( rapl.model == 45 ) {
    rapl.read_pp1 = 0;
    rapl.read_dram = 1; /* XXX: check BIOS option. memory optimization */
  } else if( rapl.model == 42 || rapl.model == 58 ) {
    rapl.read_pp1 = 1;
    rapl.read_dram = 0;
  } else {
    printf("\n");
    printf("=== WARNING ===\n");
    printf("\n");
    printf("RAPL may not be supported in your cpu\n");
    printf("Please e-mail to kazutomo@mcs.anl.gov,\n");
    printf("including your machine info(e.g. /proc/cpuinfo)\n");
    printf("\n");
  }

  /* Calculate the units used */
  val=_read_msr(rapl.pkg[0].fd ,MSR_RAPL_POWER_UNIT);
  
  rapl.power_units  = pow(0.5,(double)( val     &0xf ));
  rapl.energy_units = pow(0.5,(double)((val>>8) &0x1f));
  rapl.time_units   = pow(0.5,(double)((val>>16)&0xf ));

  for(i=0;i<rapl.npkg;i++) {
    rapl.pkg[i].pkg_power_info  = _read_msr(rapl.pkg[i].fd,MSR_PKG_POWER_INFO );  
    if(rapl.read_dram)
      rapl.pkg[i].dram_power_info = _read_msr(rapl.pkg[i].fd,MSR_DRAM_POWER_INFO);  

    rapl.pkg[i].pkg_power_limit = _read_msr(rapl.pkg[i].fd,MSR_PKG_POWER_LIMIT);
    if(rapl.read_dram)
      rapl.pkg[i].dram_power_limit= _read_msr(rapl.pkg[i].fd,MSR_DRAM_POWER_LIMIT);
    rapl.pkg[i].pp0_power_limit = _read_msr(rapl.pkg[i].fd,MSR_PP0_POWER_LIMIT);
    if(rapl.read_pp1) 
      rapl.pkg[i].pp1_power_limit = _read_msr(rapl.pkg[i].fd,MSR_PP1_POWER_LIMIT);
#if 0
    rapl.pkg[i].pp0_policy  = _read_msr(rapl.pkg[i].fd,MSR_PP0_POLICY);
    if(rapl.read_pp1) 
      rapl.pkg[i].pp1_policy  = _read_msr(rapl.pkg[i].fd,MSR_PP1_POLICY);
#endif
  }
  return 0;
}

void rapl_reader_print_info_header(FILE* fp)
{
  int i;
  fprintf(fp,"# CPU_MODEL=%d\n", rapl.model);
  fprintf(fp,"# N_SOCKETS=%d\n", rapl.npkg);
  if( rapl.read_dram) 
    fprintf(fp,"# DRAM_DOMAIN=yes\n");
  if( rapl.read_pp1) 
    fprintf(fp,"# PP1_DOMAIN=yes\n");

  fprintf(fp, "# POWER_UNITS=%lf\n",  rapl.power_units );
  fprintf(fp, "# ENERGY_UNITS=%lf\n", rapl.energy_units );
  fprintf(fp, "# TIME_UNITS=%lf\n",   rapl.time_units   );


  for(i=0; i<rapl.npkg; i++) {
    uint64_t tmpval;
    fprintf(fp,"# SOCKET%d_COREID=%d\n", i, rapl.pkg[i].coreid);

    tmpval = rapl.pkg[i].pkg_power_info;
    fprintf(fp,"# SOCKET%d_POWER_INFO: THERMAL_SPEC=%.1f MIN=%.1f MAX=%.1f TIME=%e\n", i,  
	    _get_power_info_thermal_spec(tmpval),
	    _get_power_info_min_power(tmpval),
	    _get_power_info_max_power(tmpval),
	    _get_power_info_time_window(tmpval) );

    tmpval = rapl.pkg[i].pkg_power_limit;
    fprintf(fp,"# SOCKET%d_POWER_LIMIT: LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LOCK=%d\n",i,
	    _get_power_limit_power1(tmpval),
	    _get_power_limit_enable1(tmpval),
	    _get_power_limit_clamp1(tmpval),
	    _get_power_limit_time1(tmpval),
	    _get_power_limit_power2(tmpval),
	    _get_power_limit_enable2(tmpval),
	    _get_power_limit_clamp2(tmpval),
	    _get_power_limit_time2(tmpval),
	    _get_power_limit_lock(tmpval)
	    );
    tmpval = rapl.pkg[i].pkg_power_limit_new;
    if(tmpval>0) {
      fprintf(fp,"# SOCKET%d_POWER_LIMIT: LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LIMIT1=%.1f E1=%d C1=%d TIME1=%e  LOCK=%d\n",i,
	      _get_power_limit_power1(tmpval),
	      _get_power_limit_enable1(tmpval),
	      _get_power_limit_clamp1(tmpval),
	      _get_power_limit_time1(tmpval),
	      _get_power_limit_power2(tmpval),
	      _get_power_limit_enable2(tmpval),
	      _get_power_limit_clamp2(tmpval),
	      _get_power_limit_time2(tmpval),
	      _get_power_limit_lock(tmpval)
	      );
    }
  }
}

void rapl_reader_print_info_footer(FILE* fp)
{
  int i;
  double total_sampling=0.0;
  
  
  for(i=0; i<rapl.npkg; i++) {
    int lastpos = 0;
    if( rapl.pkg[i].pos == 0 ) {
      lastpos = 1;
    }
    fprintf(fp,"# SOCKET%d_ELAPSED=%f\n", i, rapl.pkg[i].timestamp[lastpos]-rapl.pkg[i].timestamp_start );

    fprintf(fp,"# SOCKET%d_PKG_ENERGY=%lf\n", i, rapl.energy_units*rapl.pkg[i].e_total.pkg );
    fprintf(fp,"# SOCKET%d_PP0_ENERGY=%lf\n", i, rapl.energy_units*rapl.pkg[i].e_total.pp0 );
    total_sampling += rapl.energy_units*rapl.pkg[i].e_total.pkg;

    if( rapl.read_pp1) {
      fprintf(fp,"# SOCKET%d_PP1_ENERGY=%lf\n", i, rapl.energy_units*rapl.pkg[i].e_total.pp1 );
      total_sampling += rapl.energy_units*rapl.pkg[i].e_total.pp1;
    }
    if( rapl.read_dram) {
      fprintf(fp,"# SOCKET%d_DRAM_ENERGY=%lf\n", i, rapl.energy_units*rapl.pkg[i].e_total.mem );
      total_sampling += rapl.energy_units*rapl.pkg[i].e_total.mem;
    }
  }
  fprintf(fp,"# TOTAL_PKG_ENERGY=%lf\n",   total_sampling);
}



void rapl_reader_finalize(void)
{
  int i;

  rapl_reset_power_limit(); /* back to original value */
  /* we probably need to trap interrupt signal, etc */

  for(i=0;i<rapl.npkg;i++) {
    close(rapl.pkg[i].coreid);
  }    
}


static  uint32_t _read_energy(int pkgid, int offset)
{
  uint32_t v;
  v = (uint32_t)(_read_msr(rapl.pkg[pkgid].fd,offset)&0xffffffff);
  return v;
}
#if 0
static double _read_energy(int pkgid, int offset)
{
  double v=0.0;
  int i;
  for(i=0;i<10;i++) { /* retry */
    v = (double)_read_msr(rapl.pkg[pkgid].fd,offset)*rapl.energy_units;
    if(v>=0.0) break;
    usleep(1);
  }
  if(v<0.0) {
    printf("Warning: data corruption!\n");
    v = 0.0;
  }
  return v;
}
#endif

void rapl_reader_snap(void) 
{
  int i;
  /* initial energy snapshot */
  for(i=0;i<rapl.npkg;i++) {
    rapl_pkg_t *p = &(rapl.pkg[i]);

    p->timestamp[p->pos] = gettimesec();
    //printf("ts%d %lf\n",p->pos,  p->timestamp[p->pos]);
    p->e[p->pos].pkg = _read_energy(i,MSR_PKG_ENERGY_STATUS);  
    p->e[p->pos].pp0 = _read_energy(i,MSR_PP0_ENERGY_STATUS);  
    if(rapl.read_pp1) 
      p->e[p->pos].pp1 = _read_energy(i,MSR_PP1_ENERGY_STATUS);  
    if(rapl.read_dram)
      p->e[p->pos].mem = _read_energy(i,MSR_DRAM_ENERGY_STATUS);

    if(p->count==0) {
      p->e_start = p->e[p->pos];
      p->timestamp_start = p->timestamp[p->pos];
      p->e_total.pkg = p->e_total.pp0 = p->e_total.pp1 = p->e_total.mem = 0.0;
    } else {
      int lastpos = 0;;
      if(p->pos==0) lastpos = 1;
      p->e_total.pkg += (p->e[p->pos].pkg - p->e[lastpos].pkg);
      p->e_total.pp0 += (p->e[p->pos].pp0 - p->e[lastpos].pp0);
      p->e_total.pp1 += (p->e[p->pos].pp1 - p->e[lastpos].pp1);
      p->e_total.mem += (p->e[p->pos].mem - p->e[lastpos].mem);
    }

    //    p->pos++;
    //    if(p->pos > 1) p->pos = 0;
    if (p->pos)
	    p->pos = 0;
    else
	    p->pos = 1;

    p->count ++;

  }
}


void rapl_reader_start(void) 
{
  int i;
  for(i=0;i<rapl.npkg;i++) {
    rapl.pkg[i].count = 0;
  }
  rapl_reader_snap();
}

void rapl_reader_stop(void) 
{
  rapl_reader_snap();
}

void rapl_reader_report_delta(void) 
{
  int i;
  int ep,sp;
  double total=0.0;
  double subtotal = 0.0;

  if( rapl.pkg[0].count < 2 ) {
    printf("warning: no measured data\n");
    return ;
  }

  if( rapl.pkg[0].pos == 1 ) {
    ep = 0;
    sp = 1;
  } else {
    ep = 1;
    sp = 0;
  }
    
  
  //  printf("rapl.start=%f sec\n",rapl.pkg[0].timestamp_start);
  printf("rapl.elapsed=%f sec\n",rapl.pkg[0].timestamp[ep]-rapl.pkg[0].timestamp[sp]);

  /*
      p->e_total.pkg
      p->e_total.pp0
      p->e_total.pp1
      p->e_total.mem
  */
  
  for(i=0;i<rapl.npkg;i++) {
	  subtotal = 0.0;
    //    printf("rapl.pkg%d.pkg=%g J\n", i, rapl.pkg[i].e[ep].pkg-rapl.pkg[i].e[sp].pkg);
    printf("rapl.pkg%d.pkg=%g J\n", i, rapl.energy_units*rapl.pkg[i].e_total.pkg);
    printf("rapl.pkg%d.pp0=%g J\n", i, rapl.energy_units*rapl.pkg[i].e_total.pp0);
    subtotal +=  rapl.energy_units*rapl.pkg[i].e_total.pkg;

    if(rapl.read_pp1)  {
      printf("rapl.pkg%d.pp1=%g J\n", i, rapl.energy_units*rapl.pkg[i].e_total.pp1);
      subtotal +=   rapl.energy_units*rapl.pkg[i].e_total.pp1;
    }

    if(rapl.read_dram) {
      printf("rapl.pkg%d.mem=%g J\n", i, rapl.energy_units*rapl.pkg[i].e_total.mem);
      subtotal += rapl.energy_units*rapl.pkg[i].e_total.mem;
    }
    total += subtotal;
    printf("rapl.pkg%d.total=%g J\n", i, subtotal);
  }
  printf("rapl.total=%g J\n", total);

}

const rapl_reader_t*  rapl_reader_get_struct(void)
{
  return (const rapl_reader_t*)&rapl;
}


int rapl_reader_get_power(double *timestamp, double *pkg, double *pp0, 
			 double *pp1, double *mem )
{
  int i;
  int ep,sp;
  double delta;


  if( rapl.pkg[0].pos == 1 ) {
    ep = 0;
    sp = 1;
  } else {
    ep = 1;
    sp = 0;
  }

  if(!timestamp) return 0;

  *timestamp=rapl.pkg[0].timestamp[ep];

  if( !pkg ) {
    return rapl.npkg;
  }

  delta = rapl.pkg[0].timestamp[ep]-rapl.pkg[0].timestamp[sp];


  if(pkg) {
    for(i=0;i<rapl.npkg;i++) {
      pkg[i] = (rapl.pkg[i].e[ep].pkg-rapl.pkg[i].e[sp].pkg)/delta;
      pkg[i] *= rapl.energy_units;
    }
  }
  if(pp0) {
    for(i=0;i<rapl.npkg;i++) {
      pp0[i] = (rapl.pkg[i].e[ep].pp0-rapl.pkg[i].e[sp].pp0)/delta;
      pp0[i] *= rapl.energy_units;
    }
  }
  if(pp1) {
    if(rapl.read_pp1)  {
      for(i=0;i<rapl.npkg;i++) {
	pp1[i] = (rapl.pkg[i].e[ep].pp1-rapl.pkg[i].e[sp].pp1)/delta;
	pp1[i] *= rapl.energy_units;
      }
    } else {
      pp1[0] = -1.0;
    }
  }
  if(mem) {
    if(rapl.read_dram)  {
      for(i=0;i<rapl.npkg;i++) {
	mem[i] = (rapl.pkg[i].e[ep].mem-rapl.pkg[i].e[sp].mem)/delta;
	mem[i] *= rapl.energy_units;
      }
    } else {
      mem[0] = -1.0;
    }
  }
  return rapl.npkg;
}



int rapl_reader_get_energy(double *_delta, double *pkg, double *pp0, 
			 double *pp1, double *mem )
{
  int i;
  int ep,sp;
  double delta;


  if( rapl.pkg[0].pos == 1 ) {
    ep = 0;
    sp = 1;
  } else {
    ep = 1;
    sp = 0;
  }


  if( !pkg ) {
    return rapl.npkg;
  }

  delta = rapl.pkg[0].timestamp[ep]-rapl.pkg[0].timestamp[sp];
  *_delta = delta;

  if(pkg) {
    for(i=0;i<rapl.npkg;i++) {
	    pkg[i] = (rapl.pkg[i].e[ep].pkg-rapl.pkg[i].e[sp].pkg);
      pkg[i] *= rapl.energy_units;
    }
  }
  if(pp0) {
    for(i=0;i<rapl.npkg;i++) {
	    pp0[i] = (rapl.pkg[i].e[ep].pp0-rapl.pkg[i].e[sp].pp0);
      pp0[i] *= rapl.energy_units;
    }
  }
  if(pp1) {
    if(rapl.read_pp1)  {
      for(i=0;i<rapl.npkg;i++) {
	      pp1[i] = (rapl.pkg[i].e[ep].pp1-rapl.pkg[i].e[sp].pp1);
	pp1[i] *= rapl.energy_units;
      }
    } else {
      pp1[0] = -1.0;
    }
  }
  if(mem) {
    if(rapl.read_dram)  {
      for(i=0;i<rapl.npkg;i++) {
	      mem[i] = (rapl.pkg[i].e[ep].mem-rapl.pkg[i].e[sp].mem);
	mem[i] *= rapl.energy_units;
      }
    } else {
      mem[0] = -1.0;
    }
  }
  return rapl.npkg;
}


/* for fortran binding */
void  rapl_reader_init_(void)
{
  rapl_reader_init();
}
void  rapl_reader_snap_(void)
{
  rapl_reader_snap();
}

void  rapl_reader_start_(void)
{
  rapl_reader_start();
}
void  rapl_reader_stop_(void)
{
  rapl_reader_stop();
}
void  rapl_reader_report_delta_(void)
{
  rapl_reader_report_delta();
}

void  rapl_reader_finalize_(void)
{
  rapl_reader_finalize();
}

