/* for timing 10000 sweeps of n=20 at beta=6.

compile: 
gcc -O2 HM_time.c -lm -o HM_time

run:
./HM_time

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NMAX 100 // maximum matrix size
#define ORTHIZE 1000 // re-orthogonalize after this number of sweeps

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
  
  
void orthogonalize()
  {
  int i1,i2,j;
  double h11,h12;

  for(i1=0;i1<n;++i1)
    {
    for(j=0;j<n;++j)
      u[j]=h[i1][j];
		
    for(i2=0;i2<i1;++i2)
      {
      h12=.0;
      for(j=0;j<n;++j)
        h12+=u[j]*h[i2][j];
			
      h12/=n;
		
      for(j=0;j<n;++j)
        u[j]-=h12*h[i2][j];
      }
		
    h11=.0;
    for(j=0;j<n;++j)
      h11+=u[j]*u[j];
		
    h11=sqrt(h11/n);
    for(j=0;j<n;++j)
      h[i1][j]=u[j]/h11;
    }
  }
  
  
void setenergy()
{
int i,j;

energy=.0;

for(i=0;i<n;++i)
for(j=0;j<n;++j)
  energy+=-fabs(h[i][j]);
}

  
double sweep()
  {
  int i,j;
  double a;
  static int sweepcount=0;
  
  if(++sweepcount==ORTHIZE)
    {
    orthogonalize();
    setenergy();
    
    sweepcount=0;
    }
    
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
    
  setenergy();
  	
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
  

int main(int argc,char* argv[])
  {
  int sweeps;
  
  n=20;
  beta=6.;
  sweeps=10000;
  
  srand(time(NULL));
  
  init(sweeps);
  
  return 0;
  }
  
  
  
  