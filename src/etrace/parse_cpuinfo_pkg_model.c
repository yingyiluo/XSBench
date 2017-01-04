/*
  find a list of uniq core per pkg

  Kazutomo Yoshii <kazutomo@mcs.anl.gov>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

int parse_cpuinfo_pkg_model(int *cores, int ncores,int *model)
{
  FILE* fp;
  char buf[128];
  regex_t r1, r2, r3;
  regmatch_t pm[3];
  char *s1 = "processor[^[:digit:]]+([[:digit:]]+)";
  char *s2 = "physical id[^[:digit:]]+([[:digit:]]+)";
  char *s3 = "model[[:blank:]]+:[[:blank:]]+([[:digit:]]+)";
  int i, rc;
  int coreid=-1, pkgid=-1;
  char tmpbuf[10];
  int tmplen;
  int pkgid_max=-1;

  *model = 0;

  for(i=0;i<ncores;i++)  cores[i] = -1;


  rc = regcomp(&r1, s1, REG_EXTENDED);
  if( rc ) {
    fprintf(stderr, "failed to compile %s\n",s1);   exit(1);
  }
  rc = regcomp(&r2, s2, REG_EXTENDED);
  if( rc ) {
    fprintf(stderr, "failed to compile %s\n",s2);   exit(1);
  }
  rc = regcomp(&r3, s3, REG_EXTENDED);
  if( rc ) {
    fprintf(stderr, "failed to compile %s\n",s3);   exit(1);
  }


  fp = fopen( "/proc/cpuinfo", "r" );

  while( fgets(buf, sizeof(buf), fp ) != (char*)0 ) {

    rc = regexec(&r1, buf, 3,  pm, 0);
    if( rc == 0 ) {
      if( pm[1].rm_so == -1 ) {
	printf("failed to match: %s\n", s1);
	exit(1);
      }
      tmplen = pm[1].rm_eo-pm[1].rm_so+1;
      if(tmplen>5) {
	printf("processor id is too long %s:\n", buf);
	exit(1);
      }
      strncpy(tmpbuf,buf+pm[1].rm_so, tmplen);
      tmpbuf[tmplen] = 0;
      coreid = atoi(tmpbuf);
    }

    rc = regexec(&r2, buf, 3,  pm, 0);
    if( rc == 0 ) {
      if( pm[1].rm_so == -1 ) {
	printf("failed to match: %s\n", s2);
	exit(1);
      }
      tmplen = pm[1].rm_eo-pm[1].rm_so+1;
      if(tmplen>2) {
	printf("physical id is too long: %s", buf);
	exit(1);
      }
      strncpy(tmpbuf,buf+pm[1].rm_so, tmplen);
      tmpbuf[tmplen] = 0;
      pkgid = atoi(tmpbuf);
      
      /* XXX: assume physical id is continuous .
	 "processor id" is parsed before "physical id".
	 Use the first appeared core id as pkgcore.
       */
      if(pkgid>=ncores) {
	printf("physical id is too big: %d\n",pkgid);
	exit(1);
      }
      if( cores[pkgid]==-1 ) {

	if( coreid==-1 ) {
	  printf("coreid is not parsed yet\n");
	  exit(1);
	}

	cores[pkgid] = coreid;
	if( pkgid> pkgid_max ) pkgid_max = pkgid;
      }
    }

    if( *model == 0 ) {
      rc = regexec(&r3, buf, 3,  pm, 0);
      if( rc == 0 ) {
	if( pm[1].rm_so == -1 ) {
	  printf("failed to match: %s\n", s2);
	  exit(1);
	}

	tmplen = pm[1].rm_eo-pm[1].rm_so+1;
	if(tmplen>3) {
	  printf("model is too long: %s", buf);
	  exit(1);
	}
	strncpy(tmpbuf,buf+pm[1].rm_so, tmplen);
	*model = atoi(tmpbuf);
      }
    }
  }
  fclose(fp);

  return pkgid_max+1;
}

#ifdef __TEST_MAIN__
int main()
{
  int pkgcores[8];
  int npkgcores;
  int i;
  int model;
  npkgcores = parse_cpuinfo_pkg_model(pkgcores,8,&model);

  printf("model=%d\n", model);
  for(i=0;i<npkgcores;i++) {
    printf("%d: %d\n", i, pkgcores[i]);
  }
  return 0;
}
#endif
