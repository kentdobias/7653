
#pragma once
#include "randutils/randutils.hpp"
#include <vector>

inline double Ei(double x) { return -fabs(x); }

class Orthogonal {
private:
  unsigned d;
  std::vector<double> m;

public:
  Orthogonal(unsigned size) : m(size * size) {
    d = size;
    for (unsigned i = 0; i < size; i++) {
      m[size * i + i] = sqrt(size);
    }
  }

  unsigned size() const { return d; }

  double& operator()(unsigned i, unsigned j) { return m[d * i + j]; }

  const double& operator()(unsigned i, unsigned j) const { return m[d * i + j]; }

  double energy() const {
    double E = 0;

    for (unsigned i = 0; i < this->size(); i++) {
      for (unsigned j = 0; j < this->size(); j++) {
        E += Ei(this->operator()(i, j));
      }
    }

    return E;
  }
};

class Givens {
private:
  Orthogonal& m;

  bool transpose;
  unsigned axis_1;
  unsigned axis_2;
  double Δθ;

  std::vector<double> rows;

public:
  Givens(Orthogonal& m, bool t, unsigned a1, unsigned a2, double θ0, randutils::mt19937_rng& rng)
      : m(m), rows(2 * m.size()) {
    transpose = t;
    axis_1 = a1;
    axis_2 = a2;
    Δθ = rng.uniform(-θ0, θ0);
  }

  Givens(Orthogonal& m, double θ0, randutils::mt19937_rng& rng) : m(m), rows(m.size()) {
    Δθ = rng.uniform(-θ0, θ0);
    unsigned axis1axis2 = rng.uniform((unsigned)0, m.size() * (m.size() - 1) - 1);

    axis_1 = axis1axis2 / (m.size() - 1);
    axis_2 = axis1axis2 % (m.size() - 1);
    transpose = axis_2 >= axis_1;
    if (transpose) {
      axis_2++;
    }
  }

  double tryRotation() {
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

      ΔE -= Ei(m1i) + Ei(m2i);

      m1i_new = c * m1i + s * m2i;
      m2i_new = c * m2i - s * m1i;

      ΔE += Ei(m1i_new) + Ei(m2i_new);

      rows[i] = m1i_new;
      rows[m.size() + i] = m2i_new;
    }

    return ΔE;
  }

  void acceptRotation() const {
    for (unsigned i = 0; i < m.size(); i++) {
      if (transpose) {
        m(i, axis_1) = rows[i];
        m(i, axis_2) = rows[m.size() + i];
      } else {
        m(axis_1, i) = rows[i];
        m(axis_2, i) = rows[m.size() + i];
      }
    }
  }
};

class MCMC {
private:
  randutils::mt19937_rng rng;

public:
  double θ0;
  double β;
  double E;
  Orthogonal M;

  MCMC(unsigned n, double β0) : M(n) {
    β = β0;
    θ0 = M_PI;
    E = M.energy();
  }

  bool step(Givens& g) {
    double ΔE = g.tryRotation();

    if (ΔE < 0 || exp(-β * ΔE) > rng.uniform((double)0.0, 1.0)) {
      E += ΔE;
      g.acceptRotation();
      return true;
    } else {
      return false;
    }
  }

  void tune(unsigned N, double ε) {
    for (unsigned i = 0; i < N; i++) {
      Givens g(M, θ0, rng);
      bool stepAccepted = this->step(g);
      if (stepAccepted) {
        θ0 *= 1 + ε;
      } else {
        θ0 /= 1 + ε;
      }
    }
  }

  void sweep() {
    for (unsigned i = 0; i < M.size() - 1; i++) {
      for (unsigned j = i + 1; j < M.size(); j++) {
        Givens g1(M, false, i, j, θ0, rng);
        this->step(g1);
        Givens g2(M, true, i, j, θ0, rng);
        this->step(g2);
      }
    }
  }

  void run(unsigned N) {
    for (unsigned i = 0; i < N; i++) {
      this->sweep();
    }
  }
};

