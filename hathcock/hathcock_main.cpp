#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <typeinfo>

// Achknowledgements: Credit to Alen Senanian for tips regarding c++ optimization

const double PI = acos(-1);

double rotate_rows(double ** U, double const& beta, int const& N,double const& MAX_ROT = PI/4);
double hamiltonian (double ** U, int const& N);

template <typename U>
constexpr inline U d_abs(U u)
{
    return u > 0? u : -u;
}

int main(int argc, char ** argv)
{
    constexpr int N = 12;
    constexpr int N_STEPS = 1000000;
    constexpr double beta = 5.5;
    
    srand(time(NULL));
    
    // allocate memory for orthogonal matrix
    double ** U = (double**) calloc(N,sizeof(double*));
    for (int i = 0 ; i < N; i++) U[i] = (double*) calloc(N,sizeof(double));
    
    // initialize entries
    for (int i = 0; i < N; i++) U[i][i] = 1.;
    
    double energy = hamiltonian(U,N);
    printf("Initial Energy = %f\n",energy);
    
    energy = 0;
    double energy_sq = 0;
    
    double theta = 0;
    //
    int n_meas = 0;
    
    // perform Hastings-Metropolis
    for (int i = 0; i < N_STEPS; i++)
    {
        theta = rotate_rows(U,beta,N);
        // measure energy periodically
        if (i % 100 == 0)
        {
            double energy_measure = 0;
            n_meas++;
            if (n_meas > 10)
                energy_measure = hamiltonian(U,N);
                energy += energy_measure;
                energy_sq += pow(energy_measure,2);
            
        }
    }
    
    printf("Average Energy = %f\n", energy/(n_meas-10)/N/N);
    printf("Heat Capacity = %f\n", pow(beta, 2) * (energy_sq/(n_meas-10) - pow((energy/(n_meas-10)),2)) / N / N);

    for (int i = 0; i < N; i++) free(U[i]);
    free(U);
    
    return 0;
}

double rotate_rows(double ** U, double const& beta, int const& N,
                double const& MAX_ROT)
{
    
    const int i = rand() % N;
    int j = rand() % N;
    while(i==j) j = rand() % N;
    
    const double theta = MAX_ROT * (double) (rand() /(RAND_MAX + 1.));

    const double c = cos(theta);
    const double s = sin(theta);
    
    double * U_new_i = (double*) malloc(sizeof(double)*N);
    double * U_new_j = (double*) malloc(sizeof(double)*N);
    
    /* Rotate chosen rows i and j */
    for (int k = 0; k < N; k++)
    {
        U_new_i[k] = U[i][k] * c - U[j][k] * s;
        U_new_j[k] = U[j][k] * c + U[i][k] * s;
    }
    
    double deltaH = 0;

    /*Only energy change occurs on rotated rows, add up these contributions*/
    for (int k = 0; k < N; k++)
    {
        deltaH += d_abs(U_new_i[k]) + d_abs(U_new_j[k]) - d_abs(U[i][k]) -  d_abs(U[j][k]);
    }
    deltaH *= -sqrt(N);
    
    double r = (double) (rand() /( RAND_MAX + 1.));
    bool accepted = false;
    
    if (r < exp(-beta*deltaH))
    {
        for (int k = 0; k < N; ++k)
        {
            U[i][k] = U_new_i[k];
            U[j][k] = U_new_j[k];
        }
        accepted = true;
    }
    
    free(U_new_i);
    free(U_new_j);
    
    return theta;
}

double hamiltonian (double ** U, int const& N)
{
    double accum_energy = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; ++j)
        {
            accum_energy += d_abs(U[i][j]);
        }
    }
    return -sqrt(N)*accum_energy;
}


