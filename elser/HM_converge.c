/* Hadamard-model C program to study convergence

compile: 
gcc -O2 HM_converge.c -lm -o HM_converge

run:
./HM_converge n beta outfile &

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NMAX 100
#define B 20

int n;
double h[NMAX][NMAX],u[NMAX],v[NMAX];
double energy,angle,beta,twopi,accept;


double urand()
  {
  return ((double)rand())/RAND_MAX;
  }


void rowrot(int i,int j,double a)
  {
  int k;
  double eold,enew,c,s;
  
  eold=0.;
  for(k=0;k<n;++k)
    {
    u[k]=h[i][k];
    v[k]=h[j][k];
		
    eold+=-fabs(u[k])-fabs(v[k]);
    }
	
  c=cos(a);
  s=sin(a);
  
  enew=0.;
  for(k=0;k<n;++k)
    {
    h[i][k]=c*u[k]+s*v[k];
    h[j][k]=-s*u[k]+c*v[k];
		
    enew+=-fabs(h[i][k])-fabs(h[j][k]);
    }
	
  if(exp(-beta*(enew-eold))>urand())
    {
    energy+=enew-eold;
    accept++;
    }
  else
    for(k=0;k<n;++k)
      {
      h[i][k]=u[k];
      h[j][k]=v[k];
      }
  }
  
  
void colrot(int i,int j,double a)
  {
  int k;
  double eold,enew,c,s;
  
  eold=0.;
  for(k=0;k<n;++k)
    {
    u[k]=h[k][i];
    v[k]=h[k][j];
		
    eold+=-fabs(u[k])-fabs(v[k]);
    }
	
  c=cos(a);
  s=sin(a);
  
  enew=0.;
  for(k=0;k<n;++k)
    {
    h[k][i]=c*u[k]+s*v[k];
    h[k][j]=-s*u[k]+c*v[k];
		
    enew+=-fabs(h[k][i])-fabs(h[k][j]);
    }
	
  if(exp(-beta*(enew-eold))>urand())
    {
    energy+=enew-eold;
    accept++;
    }
  else
    for(k=0;k<n;++k)
      {
      h[k][i]=u[k];
      h[k][j]=v[k];
      }
  }
  
  
double sweep()
  {
  int i,j;
  double a;
  
  accept=0;
  
  for(i=0;i<n-1;++i)
  for(j=i+1;j<n;++j)
    {
    a=angle*(urand()-.5);
    rowrot(i,j,a);
    
    a=angle*(urand()-.5);
    colrot(i,j,a);
    }
    
  return ((double)accept)/(n*(n-1));
  }
  
  
void init(int sweeps)
  {
  int i,j,s;

  for(i=0;i<n;++i)
  for(j=0;j<n;++j)
    h[i][j]=0.;
	
  for(i=0;i<n;++i)
    h[i][i]=sqrt((double)n);
    
  energy=-sqrt((double)n)*n;
  	
  twopi=4.*acos(0.);
  angle=twopi;
  
  for(s=0;s<sweeps;++s)
	{
    if(sweep()<.5)
      angle*=.99;
    else
      {
      angle*=1.01;
      if(angle>twopi)
        angle=twopi;
      }
    }
  }
  
  
void average(int sweeps,double *aenergy,double *heatcap)
  {
  int s;
  double e,e2;
  
  init(sweeps);
  
  e=0.;
  e2=0.;
  
  for(s=0;s<sweeps;++s)
    {
    sweep();
    
    e+=energy;
    e2+=energy*energy;
    }
    
  e/=sweeps;
  e2/=sweeps;
  
  *aenergy=e/(n*n);
  *heatcap=beta*beta*(e2-e*e)/(n*n);
  }
  
  
void measure(int sweeps,double *eave,double *eerr,double *cave,double *cerr)
  {
  int b;
  double e,c,eave2,cave2;
  
  *eave=0.;
  eave2=0.;
  
  *cave=0.;
  cave2=0.;
  
  for(b=0;b<B;++b)
    {
    average(sweeps,&e,&c);
    
    *eave+=e;
    eave2+=e*e;
    
    *cave+=c;
    cave2+=c*c;
    }
    
  *eave/=B;
  eave2/=B;
  *eerr=sqrt((eave2-(*eave)*(*eave))/B);
  
  *cave/=B;
  cave2/=B;
  *cerr=sqrt((cave2-(*cave)*(*cave))/B);
  }
  

int main(int argc,char* argv[])
  {
  int p,sweeps;
  char *outfile;
  double eave,eerr,cave,cerr;
  FILE *fptr;
  
  if(argc==4)
    {
    n=atoi(argv[1]);
    beta=atof(argv[2]);
    outfile=argv[3];
    }
  else
    {
    printf("expected three arguments: n, beta, outfile\n");
    return 1;
    }
    
  fptr=fopen(outfile,"w");
  fprintf(fptr,"      sweeps  mean energy                     heat capacity\n\n");
  fclose(fptr);
  
  srand(time(NULL));
  
  sweeps=16;
  
  for(p=0;p<8;++p)
    {
    measure(sweeps,&eave,&eerr,&cave,&cerr);
    
    fptr=fopen(outfile,"a");
    fprintf(fptr,"%12d%14.9lf +/-%12.9lf %14.9lf +/-%12.9lf\n",sweeps,eave,eerr,cave,cerr);
    fclose(fptr);
    
    sweeps*=4;
    }
  
  return 0;
  }
  
  
  
  