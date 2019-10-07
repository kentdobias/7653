#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

/**** this function is not in working condition as of now****/
/*function to calculate the frobenius absolute norm(FAN)*/
/* float FAN(float *M, int n){
    float sum=0;
    int i,j;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            sum = sum + fabsf(M[i][j]);
        }    
    }
    return sum;
} */

int main()
{
  int n, i, j, k, c, d, p,a=1000,q,r,N,T,count=0;
  double *arr;
  double theta,sum=0,mina,maxa,sum1, EE, EE2, C, ee, CV;
  arr=(double *) malloc(100*sizeof(int));
  double beta,Ef,Ei,Q,alpha,z;
  double iom[12][12], uom[12][12]; /* om = orthogonal matrix, u is updated, i is initial*/

/*   printf("Enter order of the matrix:\n");
  scanf("%d", &n); */

    n = 12;         /*n is 12 in the question*/
    beta = 5.5;     /*beta is 5.5 in the question*/
    N=200;
    T=10000;
    srand(time(0)); /*random seeding*/

/*toggle if you wanna give a non-identity matrix as your initial orthogonal matrix*/
/*   printf("Enter the initial orthogonal matrix:\n"); 
   
  for (c = 0; c < n; c++)
    for (d = 0; d < n; d++) 
      scanf("%f%f", &iom[c][d]); */

            for (c = 0; c < n; c++){ 
                for (d = 0; d < n; d++){
                    if (c==d){
                        iom[c][c]=1;
                        uom[c][c]=1;
                    }
                    else{
                        iom[c][d]=0;
                        uom[c][d]=0;
                    }            
                }
            }

/*I plan to use the identity matrix as my initial matrix*/
    for (q=0;q<N;q++){
        for(r=0;r<T;r++){
        /*parameters of the givens rotations are chosen randomly at all the steps*/       
                    mina=-0.011;
                    maxa=0.011;
                    i = (rand() % n);
                    j = (rand() % n);
                    if (i==j){
                        continue;
                    }
                    z=((double) (rand() % a))/((double) a);
                    theta = mina+(maxa-mina)*(2*3.14159)*z;

                /*Matrix multiplication, i.e., updated orthogonal matrix*/
                /*Only ith and jth row will change*/
            
                    sum=0;
                    for (c = 0; c < n; c++){
                        for (d = 0; d < n; d++){
                            sum = sum + fabs(iom[c][d]);
                        }    
                    }

                    Ei = -sqrt(n)*sum;

                    for (k=0;k<n;k++){
                        uom[i][k]=cos(theta)*iom[i][k]-sin(theta)*iom[j][k];
                        uom[j][k]=sin(theta)*iom[i][k]+cos(theta)*iom[j][k];
                    }

                    sum=0;

                    for (c = 0; c < n; c++){
                        for (d = 0; d < n; d++){
                            sum = sum + fabs(uom[c][d]);
                        }    
                    }

                    Ef = -sqrt(n)*sum;

                    Q = exp(-beta*(Ef-Ei));
                    alpha =((double) (rand() % a))/((double) a); 
                    /*Acceptance*/
                    if(alpha<Q){
                        for (c = 0; c < n; c++){
                            for (d = 0; d < n; d++){
                                iom[c][d]=uom[c][d];
                            }    
                        } 
                        count = count+1;             
                    }
                    else{
                        for (c = 0; c < n; c++){
                            for (d = 0; d < n; d++){
                                uom[c][d]=iom[c][d];
                            }    
                        }              
                    }
        }
        /*putting the energy in the stack*/
        arr[q]=Ei;
    }

    printf("\n");

    sum1 = 0; 
    sum = 0;

    for ( d = 0; d < N; d++){
        sum = sum + arr[d];
        sum1 = sum1 + arr[d]*arr[d];
    }

    EE = (double) sum/(double) N;
    EE2 = (double) sum1/(double) N;
    C = beta*beta*(EE2 - EE*EE);
    ee = (double) EE/(double) (n*n);
    CV = (double) C/(double) (n*n);

    printf("Average energy is %f\n Heat capacity is %f\n", EE,C);
    printf("Average energy per unit volume is %f\n Heat capacity per unit volume is %f\n", ee,CV);
    printf("\n");

    printf("\n");
    printf("%d\n", count);
    return 0; 
}