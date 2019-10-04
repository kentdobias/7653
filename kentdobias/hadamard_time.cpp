#include "hadamard_mcmc.hpp"

int main() {
  MCMC sim(20, 6.0);
  sim.run(1e4);

  return 0;
}
