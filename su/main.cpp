#include <cassert>
#include <cmath>
#include <cstring>
#include <random>
#include <cstdio>
#include "HadamardMCSim.hpp"

// params (overridable via g++ -D DIM=16)
#ifndef DIM
#define DIM 12
#endif
#ifndef BETA
#define BETA 5.5
#endif
#ifndef NEPOCHS
#define NEPOCHS 400
#endif
#ifndef NUMSWEEPS
#define NUMSWEEPS 4096
#endif
#ifndef INITEPOCHS
#define INITEPOCHS 0
#endif
// also can use #define NDEBUG, #define NOSTDEVS

typedef double doub;

// use typed consts to catch bad #define's
const int dim = DIM;
const doub beta = BETA;
const int n_epochs = NEPOCHS; // not including initial epoch
const int epoch_num_sweeps = NUMSWEEPS;
const int init_epochs = INITEPOCHS;
#ifndef NOSTDEVS
FILE *stdevs_file = fopen("stdevs.log", "w");
#endif

using namespace std;

// compute mean/variance and store in passed refs
void stats(doub *arr, int len, doub *mean, doub *vars)
{
    doub old_mean = 0;
    doub old_var = 0; // piggyback for i < 2 cases
    doub sum = 0;
    *mean = 0;
    for (int i = 0; i < len; i++)
    {
        sum += arr[i];
        old_mean = *mean;
        *mean = old_mean + (arr[i] - old_mean) / (i + 1);
        if (i < 2)
        {
            vars[i] = old_var + (arr[i] - old_mean) * (arr[i] - *mean);
            old_var = vars[i];
        }
        else
        {
            vars[i] = (vars[i - 1] + (arr[i] - old_mean) * (arr[i] - *mean))
                / (i - 1);
        }
    }
}
void run_mcmc()
{
    const doub reject_target = 0.5;

    HadamardMCSim<dim, doub> mat; // IC = identity
    doub q_max = sqrt(log(2) / (beta * sqrt(2))); // inital guess, ~5-10% error
    doub guess_q; // intermediate q_max value
    doub means[n_epochs];
    doub vars[n_epochs];
    doub mean_mean, var_mean;
    doub mean_vars[n_epochs], var_vars[n_epochs];

    // special initialization epoch, do manually while adjusting q_max
    printf("q_max init %f\n", q_max);
    for (int i = 0; i < init_epochs; i++)
    {
        mat.reset_counts();
        for (int j = 0; j < epoch_num_sweeps; j++)
        {
            mat.sweep(beta, q_max);
        }
        // if rejects / steps > reject_target, q_max decreases
        // change by smaller amount in later subepochs
        guess_q = q_max / (
                ((doub) mat.get_rejects() / mat.get_steps()) / reject_target);
        q_max = (q_max + guess_q) / 2;
    }
    printf("q_max final %f, last init epoch rejects/steps = %d/%d (%.2f%%)\n",
            q_max, mat.get_rejects(), mat.get_steps(),
            ((doub) mat.get_rejects() / mat.get_steps()) * 100);

    // do all sweeps now
    for (int i = 0; i < n_epochs; i++)
    {
        mat.reset_counts();
        for (int j = 0; j < epoch_num_sweeps; j++)
        {
            mat.sweep(beta, q_max);
        }
#ifndef NDEBUG
        printf("H = %f, (rejects/steps) = (%d/%d), e/c = (%.5f/%.5f)\n",
                mat.H(), mat.get_rejects(), mat.get_steps(),
                mat.get_mean() / pow(dim, 2),
                mat.get_var() * pow(beta, 2) / pow(dim, 2));
#endif
        means[i] = mat.get_mean();
        vars[i] = mat.get_var();
    }

    // populate stats vars
    stats(means, n_epochs, &mean_mean, mean_vars);
    stats(vars, n_epochs, &var_mean, var_vars);
    printf("e = %f +- %f, c = %f +- %f\n",
            mean_mean / pow(dim, 2),
                sqrt(mean_vars[n_epochs - 1]) / pow(dim, 2),
            var_mean * pow(beta, 2) / pow(dim, 2),
                sqrt(var_vars[n_epochs - 1]) * pow(beta, 2) / pow(dim, 2));
    printf("Last epoch reject/total = %d/%d (%.2f%%)\n",
            mat.get_rejects(), mat.get_steps(),
            ((doub) mat.get_rejects() / mat.get_steps()) * 100);
#ifndef NOSTDEVS
    // print out errors
    for (int i = 0; i < n_epochs; i++)
    {
        fprintf(stdevs_file, "%f%s",
                sqrt(mean_vars[i] / pow(dim, 2)),
                i == n_epochs - 1 ? "\n" : ",");
    }
    for (int i = 0; i < n_epochs; i++)
    {
        fprintf(stdevs_file, "%f%s",
                sqrt(var_vars[i] * pow(beta, 2) / pow(dim, 2)),
                i == n_epochs - 1 ? "\n" : ",");
    }
#endif
}

int main(int argc, const char *argv[])
{
    run_mcmc();
    return 0;
}
