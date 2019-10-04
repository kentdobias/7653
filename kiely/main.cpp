#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>

using namespace std;

//Generates random float between 0 and 1
float unitRandom()
{
    float frac=1.0/RAND_MAX;
    return rand()*frac;
}

//Generates random float betweeen -maxVal and maxVal
float randomAngle(float maxVal)
{
    float frac=1.0/RAND_MAX;
    return rand()*frac*2*maxVal - maxVal;
}

//Generates random integer between 0 and n
int randomIndex(int n)
{
    return rand()%n;
}

//The 2D array to be rotated is passed as a pointer to a pointer, size is the length of the square matrix along one dimension,
//i and j (with i<j) are the indices (start from 0) of the rows being rotated by an angle theta.
void ElementaryTransition(float **arr, int n, int& accept, float& energy, float beta, int i, int j, float thetaRange)
{
    float u[n];
    float v[n];
    float eold = 0.;
    float theta = randomAngle(thetaRange);
    for(int k = 0; k < n; k++ )
    {
        u[k] = *(*(arr + i) + k);
        v[k] = *(*(arr + j) + k);
        eold += -fabs(u[k])-fabs(v[k]);
    }
    float enew = 0.;
    float c=cos(theta);
    float s=sin(theta);
    for(int k = 0; k < n; k++ )
    {
        *(*(arr + i) + k) = c*u[k] + s*v[k];
        *(*(arr + j) + k) = c*v[k] - s*u[k];
        enew += -fabs(*(*(arr + i) + k))-fabs(*(*(arr + j) + k));
    }
    if(exp(-sqrt(n)*beta*(enew-eold))>unitRandom())
    {
        energy += sqrt(n)*(enew-eold);
        accept++;
    }
    else
    {
        for(int k = 0; k < n; k++ )
        {
            *(*(arr + i) + k) = u[k];
            *(*(arr + j) + k) = v[k];
        }
    }
}

//Also passes the 2D array as a pointer to a pointer, along with the size of the array as an integer. Alpha is the exponent in the Hamiltonian.
float Hamiltonian(float **arr, int n, float alpha)
{
    float total;
    int i, j;
    total = 0;
    for( i = 0; i < n; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            total += pow(abs(*(*(arr + j) + i)), alpha);
        }
    }
    return -1*sqrt(n)*total;
}

//Performs Givens rotations for every pair of rows and columns
void PerformSweep(float **arr, int n, int& accept, float& energy, double& avEnergy, double& avEnergySq, float beta, float thetaRange)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<i;j++)
        {
            ElementaryTransition(arr, n, accept, energy, beta, j, i, thetaRange);
            avEnergy+=energy;
            avEnergySq+=energy*energy;
        }
    }
}

//Initializes an array, passes it to the pointer, performs numBlocks time steps, each consisting of a sweep followed by blockSize elementary steps,
//and calculates average energy and specific heat capacity.
int main()
{
    srand(static_cast<unsigned int>(std::time(nullptr)));
    float ar[12][12] = {{1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},{0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
    {0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.},
    {0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.},
    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.}};
    int n = 12;
    int accept = 0;
    int iterations = 0;
    float acceptanceProb = 0.;
    float acceptanceProbs[20];
    float energy = 0.;
    float beta = 5.5;
    float maxTheta = 3.14159/11.;
    int numBlocks = 100;
    int blockSize = 1000000;
    double tempEnergyAv = 0.;
    double energyAv[numBlocks];
    double totalEnergyAv[numBlocks-1];
    double heatCapAv[numBlocks];
    double totalHeatCapAv = 0.;
    double tempEnergyStd = 0.;
    double energyStd[numBlocks];
    double totalEnergyStd[numBlocks-1];
    double totalHeatCapStd = 0.;
    int randIndex1 = 0;
    int randIndex2 = 0;

    float **array;
    array = new float *[n];
    for(int i = 0; i < n; i++)
        array[i] = new float[n];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            array[i][j] = ar[i][j];
        }
    }
    energy = Hamiltonian(array,n,1.);

    for(int w=0;w<numBlocks;w++)
    {
        iterations = 0;
        accept = 0;
        tempEnergyAv = 0.;
        tempEnergyStd = 0.;
        PerformSweep(array,n,accept,energy,tempEnergyAv,tempEnergyStd,beta,maxTheta);
        iterations+=n*(n-1)/2;
        for(int q=0;q<blockSize;q++)
        {
            randIndex1 = randomIndex(n);
            while(randIndex1==0)
                randIndex1 = randomIndex(n);
            randIndex2 = randomIndex(randIndex1);
            ElementaryTransition(array,n,accept,energy,beta,randIndex2,randIndex1,maxTheta);
            tempEnergyAv+=energy;
            tempEnergyStd+=energy*energy;
            iterations++;
        }
        energyAv[w] = tempEnergyAv/iterations;
        heatCapAv[w] = ((tempEnergyStd/iterations)-energyAv[w]*energyAv[w]);
        energyStd[w] = sqrt(heatCapAv[w]);
        heatCapAv[w] = heatCapAv[w]*beta*beta;
        acceptanceProb=1.0*accept/iterations;
        acceptanceProbs[w]=acceptanceProb;
        if(w==0 && acceptanceProb>0.5)
        {
            maxTheta = 1.05*maxTheta;
        }
        if(w==0 && acceptanceProb<0.5)
        {
            maxTheta = 0.95*maxTheta;
        }
    }
    for(int i=0;i<numBlocks-1;i++)
    {
        if(i==0)
        {
            totalEnergyAv[i]=energyAv[i+1];
            totalEnergyStd[i]=energyStd[i+1];
        }
        else
        {
            float tempAv = 0.;
            float tempAvSq = 0.;
            float tempVarAv = 0.;
            for(int j=0;j<=i;j++)
            {
                tempAv+=energyAv[j+1]/(i+1);
                tempAvSq+=energyAv[j+1]*energyAv[j+1]/(i+1);
                tempVarAv+=energyStd[j+1]*energyStd[j+1]/(i+1);
            }
            totalEnergyAv[i]=tempAv;
            totalEnergyStd[i]=sqrt((tempVarAv)+(tempAvSq-tempAv*tempAv)/(i+1))/sqrt(i+1);
        }
    }
    for(int i=1;i<numBlocks-1;i++)
    {
        totalHeatCapAv+=heatCapAv[i+1]/(numBlocks-1);
        totalHeatCapStd+=heatCapAv[i+1]*heatCapAv[i+1]/(numBlocks-1);
    }
    totalHeatCapStd=sqrt(totalHeatCapStd-totalHeatCapAv*totalHeatCapAv)/sqrt(numBlocks-1);


//    for(int i=0;i<n-1;i++)
//    {
//        printf("%f\t", totalEnergyAv[i]);
//    }
//    printf("\n");
//    for(int i=0;i<n-1;i++)
//    {
//        printf("%f\t", heatCapAv[i]);
//    }

    cout<<"The average energy is " << totalEnergyAv[numBlocks-2] << " with a standard deviation of " << totalEnergyStd[numBlocks-2] << ".\n";
    cout<<"The average heat capacity is " << totalHeatCapAv << " with a standard deviation of " << totalHeatCapStd << ".\n";


    return 0;
}
