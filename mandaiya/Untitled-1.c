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
  int n, i, j, k, c, d, p,a=1000,q,r,N,T;
  double *arr,x;
  double theta,sum=0,mina,maxa;
  arr=(double *) malloc(100*sizeof(int));
  float beta,Ef,Ei,Q,alpha,z;
  float GRM[12][12], iom[12][12], uom[12][12]; /* om = orthogonal matrix, u is updated, i is initial*/

/*   printf("Enter order of the matrix:\n");
  scanf("%d", &n); */

    n = 12;         /*n is 12 in the question*/
    beta = 5.5;     /*beta is 5.5 in the question*/
    N=20;
    T=30;
    srand(time(0)); /*random seeding*/
    x=1000000000000000;
/*toggle if you wanna give a non-identity matrix as your initial orthogonal matrix*/
/*   printf("Enter the initial orthogonal matrix:\n"); 
   
  for (c = 0; c < n; c++)
    for (d = 0; d < n; d++) 
      scanf("%f%f", &iom[c][d]); */

            for (c = 0; c < n; c++){ 
                for (d = 0; d < n; d++){
                    if (c==d){
                        iom[c][c]=(long double) 1*x;
                        uom[c][c]=(long double) 1*x;
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
            
                    for (c = 0; c < n; c++){ 
                        for (d = 0; d < n; d++){
                            if (c==d){
                                GRM[c][c]=1;
                            }
                            else{
                                GRM[c][d]=0;
                            }            
                        }
                    }
                    mina=-0.01;
                    maxa=0.01;
                    i = (rand() % n) + 1;
                    j = (rand() % n) + 1;
                    z=((double) (rand() % a))/((double) a);
                    theta = mina+(maxa-mina)*(2*3.14159)*z;

                    GRM[i][i] = cos(theta);
                    GRM[j][j] = cos(theta);
                    GRM[i][j] = -sin(theta);
                    GRM[j][i] = sin(theta); 

                /*Matrix multiplication, i.e., updated orthogonal matrix*/
                /*Only ith and jth row will change*/
                
                    for (k=0;k<n;k++){
                        uom[i][k]=((long double) cos(theta)*(long double) iom[i][k]-(long double) sin(theta)*(long double)iom[j][k]);
                        uom[j][k]=((long double) sin(theta)*(long double) iom[i][k]+(long double) cos(theta)*(long double)iom[j][k]);
                    }

                    sum=0;
                    for (c = 0; c < n; c++){
                        for (d = 0; d < n; d++){
                            sum = (long double) sum + (long double) fabsl((long double) iom[c][d]);
                        }    
                    }

                    Ei = -sqrt(n)*sum/((double) x);

                    sum=0;

                    for (c = 0; c < n; c++){
                        for (d = 0; d < n; d++){
                            sum = (long double) sum + (long double) fabsl((long double) uom[c][d]);
                        }    
                    }

                    Ef = -sqrt(n)*sum/((double) x);    
                    Q = exp(-beta*(Ef-Ei));
                    alpha =((double) (rand() % a))/((double) a); 
                    /*Acceptance*/
                    if(alpha<Q){
                        for (c = 0; c < n; c++){
                            for (d = 0; d < n; d++){
                                iom[c][d]=uom[c][d];
                            }    
                        }              
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

    for ( d = 0; d < N; d++){
        printf("%f\t", arr[d]);
    }
    printf("\n");
    return 0; 
}