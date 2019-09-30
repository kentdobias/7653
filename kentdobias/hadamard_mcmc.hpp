/* Hadamard-model C++ Monte Carlo interface.
 *
 * The Orthogonal class contains virtual functions that must be defined in a
 * child class! Once done, link into a C++ program, construct a MCMC object
 * with matrix size and β, and use its member functions tune and run.
 */

#include "randutils/randutils.hpp"
#include <cmath>
#include <vector>

inline double Ei(double x, double α = 1) { return -pow(fabs(x), α); }

class Orthogonal {
public:
  virtual unsigned size() const = 0;
  virtual double& operator()(unsigned, unsigned) = 0;
  virtual const double& operator()(unsigned, unsigned) const = 0;

  double energy(double α = 1) const {
    double E = 0;

    for (unsigned i = 0; i < this->size(); i++) {
      for (unsigned j = 0; j < this->size(); j++) {
        E += Ei(this->operator()(i, j), α);
      }
    }

    return E;
  }
};

class Givens {
private:
  bool transpose;
  unsigned axis_1;
  unsigned axis_2;
  double Δθ;

  std::vector<double> row_1;
  std::vector<double> row_2;

public:
  Givens(unsigned size, double θ0, randutils::mt19937_rng& rng) : row_1(size), row_2(size) {
    Δθ = rng.uniform(-θ0, θ0);
    unsigned axis1axis2 = rng.uniform((unsigned)0, size * (size - 1) - 1);

    axis_1 = axis1axis2 / (size - 1);
    axis_2 = axis1axis2 % (size - 1);
    transpose = axis_2 >= axis_1;
    if (transpose) {
      axis_2++;
    }
  }

  double initializeRotation(const Orthogonal& m, double α = 1) {
    double ΔE = 0;
    double c = cos(Δθ);
    double s = sin(Δθ);

    for (unsigned i = 0; i < m.size(); i++) {
      double m1i, m2i, m1i_new, m2i_new;

      if (transpose) {
        m1i = m(i, axis_1);
        m2i = m(i, axis_2);
      } else {
        m1i = m(axis_1, i);
        m2i = m(axis_2, i);
      }

      ΔE -= Ei(m1i, α) + Ei(m2i, α);

      m1i_new = c * m1i + s * m2i;
      m2i_new = c * m2i - s * m1i;

      ΔE += Ei(m1i_new, α) + Ei(m2i_new, α);

      row_1[i] = m1i_new;
      row_2[i] = m2i_new;
    }

    return ΔE;
  }

  void acceptRotation(Orthogonal& m) const {
    for (unsigned i = 0; i < m.size(); i++) {
      if (transpose) {
        m(i, axis_1) = row_1[i];
        m(i, axis_2) = row_2[i];
      } else {
        m(axis_1, i) = row_1[i];
        m(axis_2, i) = row_2[i];
      }
    }
  }
};

template <class OrthogonalClass> class MCMC {
private:
  randutils::mt19937_rng rng;

public:
  double θ0;
  double β;
  double E;
  OrthogonalClass M;

  MCMC(unsigned n, double β0) : M(n) {
    β = β0;
    θ0 = M_PI;
    E = M.energy();
  }

  bool step() {
    Givens trial(M.size(), θ0, rng);
    double ΔE = trial.initializeRotation(M);

    if (exp(-β * ΔE) > rng.uniform((double)0.0, 1.0)) {
      E += ΔE;
      trial.acceptRotation(M);
      return true;
    } else {
      return false;
    }
  }

  void tune(unsigned N, double ε) {
    for (unsigned i = 0; i < N; i++) {
      bool stepAccepted = this->step();
      if (stepAccepted) {
        θ0 += ε;
      } else {
        θ0 -= ε;
      }
    }
  }

  void run(unsigned N) {
    for (unsigned i = 0; i < N; i++) {
      this->step();
    }
  }
};

