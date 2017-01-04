/*
  An enery trace code 

  support spawn mode(e.g. strace) and polling mode that priodically 
  print power usage.

  currently only support RAPL
  
  Kazutomo Yoshii <kazutomo@mcs.anl.gov>
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>


#include "rapl_reader.h"

#define ETRACE_VERSION "0.5"

static char *program_name;

char *_basename(char *path)
{
  char *base = strrchr(path, '/');
  return base ? base+1 : path;
}

/*
  XXX: clean up this later
*/
#define LOCKFN  "/var/tmp/etrace.lock"
static void  touch_lock(void) 
{
  int hackfd;
  hackfd = open(LOCKFN, O_WRONLY|O_CREAT|O_NOCTTY|O_NONBLOCK, 0666);
  close(hackfd);
}

static int access_lock(void)
{
  return  access(LOCKFN, R_OK);
}
static void unlink_lock(void)
{
  unlink(LOCKFN);
}
static void sigh(int signum) 
{
  unlink(LOCKFN);
  exit(1);
}



/*
  if child_pid > 0, it doesn't use timeout and waits until the child process is terminated
*/
static void  sampling_loop(FILE* fp, double interval, double timeout, int child_pid, int verbose, int relative_time)
{
  double ts, time_start, time_delta;
  double pkg[MAX_RAPL_PKG], pp0[MAX_RAPL_PKG];
  double pp1[MAX_RAPL_PKG], mem[MAX_RAPL_PKG];
  useconds_t usec;
  int status;
  int npkg;
  int i;
  int rc;
  extern double gettimesec(void);


  fprintf(fp, "# ETRACE_VERSION=%s\n",ETRACE_VERSION);

  if(verbose)
    rapl_reader_print_info_header(fp);

  rc = access_lock();
  if( rc == 0 ) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Failed to start etrace!\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Someone else is running etrace. RAPL is a sytem-wide resource.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  touch_lock();
  signal(SIGINT, sigh);


  rapl_reader_snap();
  npkg = rapl_reader_get_power(&ts,NULL,NULL, NULL,NULL);
  time_start = ts;

  if( interval>0.0 ) {
    if(relative_time)  {
      fprintf(fp, "# START_TIME=%lf\n",gettimesec());
    }

    fprintf(fp,"# since_epoch(sec)  ");
    for(i=0;i<npkg;i++) {
      fprintf(fp, "pkg%d_total(W)  pkg%d_core(W)  ", i, i);
      if( rapl_reader_domain_pp1() ) {
	fprintf(fp, "pkg%d_pp1(W)  ", i);
      }
      if( rapl_reader_domain_dram() ) {
	fprintf(fp, "pkg%d_dram(W) ", i);
      }
    }
    fprintf(fp, "\n");

    /* polling loop */
    for(;;) {
      time_delta = fmod(ts, interval);
      usec = (useconds_t)((interval-time_delta)*1000.0*1000.0);
      usleep(usec);

      rapl_reader_snap();
      rapl_reader_get_power(&ts,pkg,pp0,pp1,mem);
      if( child_pid > 0 ) {
	if( waitpid(child_pid, &status, WNOHANG) != 0 ) 
	  break;
      } else  if(ts-time_start>timeout )  {
	break;
      }

      if(relative_time) 
	fprintf(fp,"%.3lf  ", ts-time_start);
      else 
	fprintf(fp,"%.3lf  ", ts);

      for(i=0;i<npkg;i++) {
	fprintf(fp,"%lf  %lf  ", pkg[i], pp0[i]);
	if( rapl_reader_domain_pp1() ) {
	  fprintf(fp, "%lf  ", pp1[i]);
	}
	if( rapl_reader_domain_dram() ) {
	  fprintf(fp, "%lf ", mem[i]);
	}
      }
#if 0
      /* ad hoc */
      {
	extern void print_adhoc(FILE*);
	print_adhoc(fp);
      }
#endif
      fprintf(fp, "\n");
      fflush(fp);
    }
  } else {
    /* just energy usage */
    if( child_pid > 0 ) {
      waitpid(child_pid, &status, 0);
    } else {
      usleep(timeout*1000.0*1000.0);
    }
  } 

  rapl_reader_snap();

  rapl_reader_print_info_footer(fp);
}

static void usage(void)
{
  printf("\n");
  printf("Usage: %s [options] [command]\n",program_name);
  printf("\n");
  printf("[Options]\n");
  printf("\t-i sec : interval\n");
  printf("\t-t sec : timeout\n");
  printf("\t-o filename : output filename\n");
  printf("\t-r : relative time, instead of epoch\n");
  printf("\n");
  printf("\t--setplim1 watt : limit power\n");
  printf("\t--setplim2 watt : limit power (turbo)\n");
  printf("\n");
  printf("[usage examples]\n");
  printf("\n");
  printf("# spawn your program and print energy consumption \n");
  printf("$ ./%s your_program\n",program_name);
  printf("\n");
  printf("# also print power consumption every 0.5 sec\n");
  printf("$ ./%s -i 0.5  your_program\n",program_name);
  printf("\n");
  printf("# just print power consumption every 1sec for 10sec\n");
  printf("$ ./%s -i 1.0 -t 10\n",program_name);
  printf("\n");
}

int main(int argc, char *argv[])
{
  int opt;
  int verbose=0;
  int relative_time=0;
  double timeout=-1.0, interval=-1.0;
  char fn[1024];
  FILE* fp = stdout;
  double plim1=-1.0, plim2=-1.0;

  fn[0] = 0;
#if defined(USE_MSR_SAFE)
  if( access("/dev/cpu/0/msr_safe", F_OK)!=0 ) {
    printf("/dev/cpu/0/msr does not exist. Install the msr safe driver.\n");
    printf("\n");
    exit(1);
  }
#else
  if( access("/dev/cpu/0/msr", F_OK)!=0 ) {
    printf("/dev/cpu/0/msr does not exist. maybe msr is not installed?\n");
    printf("\n");
    printf("Suggestion:\n");
    printf("\n");
    printf("sudo modprobe msr\n");
    printf("\n");
    exit(1);
  }

  if( access("/dev/cpu/0/msr", R_OK)!=0 ) {
    printf("\n");
    printf("/dev/cpu/0/msr is not readable. please check the permission\n");
    printf("\n");
    printf("Suggestion:\n");
    printf("\n");
    printf("sudo chmod o+r /dev/cpu/*/msr\n");
    printf("sudo setcap cap_sys_rawio=ep etrace\n");
    printf("\n");
    exit(1);
  }
#endif

  rapl_reader_init();

  program_name = _basename(argv[0]);
  while(1) {
    int option_index = 0;
    static struct option long_options[] = {
      // name, has_arg, flag, val
      {"setplim1", 1, 0,  0 },
      {"setplim2", 1, 0,  0 },
      {0,          0, 0,  0 }
    };

    opt=getopt_long(argc, argv, "+i:t:o:rhv",
	       long_options, &option_index);
    if(opt==-1) break;
					       
    switch(opt) {
    case 0:
      switch( option_index ) {
      case 0:
	plim1 = strtod(optarg, (char**)0);
	break;
      case 1:
	plim2 = strtod(optarg, (char**)0);
	break;
      }
      break;
    case 'i':
      interval = strtod(optarg,(char**)0);
      break;
    case 't':
      timeout = strtod(optarg,(char**)0);
      break;
    case 'o':
      strncpy(fn,optarg,sizeof(fn)-1);
      fn[sizeof(fn)-1] = 0;
      break;
    case 'r':
      relative_time = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
      usage();
      exit(0);
    }
  }

  if(strlen(fn)>0) {
    fp = fopen(fn, "w");
    if(!fp) {
      perror("fopen");
      exit(1);
    }
  }

  if( plim1>0.0 && plim2>0.0 ) {
    rapl_set_power_limit(plim1, plim2);
  }

  if(argc-optind) {
    pid_t child_pid;

    /* spawn mode */
    child_pid = fork();
    if( child_pid==0 ) {
      execvp( argv[optind], argv+optind );
    } else {
      sampling_loop(fp, interval, 0, child_pid, verbose, relative_time);
    }
  } else {
    /* polling mode */
    if( timeout < 0.0 || interval < 0.0 ) {
      printf("Please specify timeout(-t) and interval(-i)\n");
      usage();
      exit(1);
    }
    sampling_loop(fp,interval, timeout, 0, verbose, relative_time);
  }

  if( fp!=stdout ) {
    fclose(fp);
  }


  rapl_reader_finalize();

  unlink_lock();


  return 0;
}


