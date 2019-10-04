/**
 * templatted matrix size so compile-time loop unroll if need be
 * - row-major order
 * - T has to be numeric (obviously)
 */
#ifndef MAT2D_HPP
#define MAT2D_HPP

// used to generate random numbers for given distributions
using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> row_or_col(0, 1);
uniform_real_distribution<> unit(0.0, 1.0);

template<int N, typename T> class HadamardMCSim
{
    private:
        T data[N * N];
        T l1_norm; // -H / sqrt(n), easier to track
        T s_n; // precompute sqrt(n)
        // temporary variables for storing candidate Givens changes before
        // acceptance/rejection step
        T cand1[N];
        T cand2[N];
        T l1_delta;

        // operator() cannot be made const, so make it private
        T& operator()(const int row, const int col)
        {
            assert(row < N);
            assert(col < N);
            return this->data[row * N + col];
        }

        /**
         * helper functions for proposing and applying Givens rotations, helpful
         * in case rotations are rejected
         *
         * do not assume cand1, cand2, l1_delta are sanitized at start of helpers
         */
        void gen_givens_row(const int row1, const int row2, const T q)
        {
            assert(row1 != row2);
            T c, s;
            s = sin(q);
            c = sqrt(1 - s * s);
            this->l1_delta = 0;
            for (int i = 0; i < N; i++)
            {
                this->cand1[i] = c * (*this)(row1, i) - s * (*this)(row2, i);
                // once assembled, this lookup is the same as storing to tmp var
                this->l1_delta += abs(this->cand1[i]);

                this->cand2[i] = s * (*this)(row1, i) + c * (*this)(row2, i);
                this->l1_delta += abs(this->cand2[i]);

                this->l1_delta -= abs((*this)(row1, i));
                this->l1_delta -= abs((*this)(row2, i));
            }
        }
        void apply_givens_row(const int row1, const int row2)
        {
            for (int i = 0; i < N; i++)
            {
                // NOTE: this gets optimized to a memcpy by g++
                (*this)(row1, i) = this->cand1[i];
                (*this)(row2, i) = this->cand2[i];
            }
            this->l1_norm += this->l1_delta;
        }
        void gen_givens_col(const int col1, const int col2, const T q)
        {
            assert(col1 != col2);
            T c, s;
            s = sin(q);
            c = sqrt(1 - s * s);
            this->l1_delta = 0;
            for (int i = 0; i < N; i++)
            {
                this->cand1[i] = c * (*this)(i, col1) + s * (*this)(i, col2);
                this->l1_delta += abs(this->cand1[i]);

                this->cand2[i] = -s * (*this)(i, col1) + c * (*this)(i, col2);
                this->l1_delta += abs(this->cand2[i]);

                this->l1_delta -= abs((*this)(i, col1));
                this->l1_delta -= abs((*this)(i, col2));
            }
        }
        void apply_givens_col(const int col1, const int col2)
        {
            for (int i = 0; i < N; i++)
            {
                (*this)(i, col1) = this->cand1[i];
                (*this)(i, col2) = this->cand2[i];
            }
            this->l1_norm += this->l1_delta;
        }

        /**********************
         * statistics section *
         **********************/
        int steps;
        int rejects;

        // track energy mean / variance online: Welford's algorithm
        T e_mean;
        T e_var_base; // = (steps - 1) * variance

        uniform_int_distribution<> get_idx;
        uniform_int_distribution<> get_idx2;

    public:
        HadamardMCSim () :
            s_n(sqrt(N)),
            steps(0),
            rejects(0),
            e_mean(0),
            e_var_base(0),
            get_idx(0, N - 1),
            get_idx2(0, N - 2)
        {
            this->l1_norm = N;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    (*this)(i, j) = i == j;
                }
            }
        }

        HadamardMCSim(const T data[N * N]) :
            s_n(sqrt(N)),
            steps(0),
            rejects(0),
            e_mean(0),
            e_var_base(0),
            get_idx(0, N - 1),
            get_idx2(0, N - 2)
        {
            this->l1_norm = 0;
            for (int i = 0; i < N * N; i++)
            {
                this->data[i] = data[i];
                this->l1_norm += abs(data[i]);
            }
        }
        ~HadamardMCSim () {}

        // accessors
        const T inline H() const { return -s_n * this->l1_norm; }
        const int get_steps() const { return this->steps; }
        const int get_rejects() const { return this->rejects; }
        const T inline get_mean() const { return this->e_mean; }
        const T inline get_var()  const
        {
            return this->e_var_base / (this->steps - 1);
        }
        // cannot overload operator() since returned T (double) is never a const
        const T get(const int row, const int col) const
        {
            assert(row < N);
            assert(col < N);
            return this->data[row * N + col];
        }

        /**
         * Proposes Givens rotations between 2 rows by angle q. Stores temporary
         * row/col/H in cand1, cand2, l1cand
         */
        void do_givens_row(const int row1, const int row2, const T q)
        {
            gen_givens_row(row1, row2, q);
            apply_givens_row(row1, row2);
        }
        /**
         * Mostly duplicated off do_givens_row() but should be faster than
         * fiddling around w/ inline helper funcs
         */
        void do_givens_col(const int col1, const int col2, const T q)
        {
            gen_givens_col(col1, col2, q);
            apply_givens_col(col1, col2);
        }

        void print()
        {
            printf("H=%f, Mat:\n", this->H());
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    printf("%f%s", (*this)(i, j), j == N - 1 ? "" : ", ");
                }
                printf("\n");
            }
        }

        /**
         * Executes an MCMC step by:
         * - generate a random Givens rotation
         * - computes whether to accept and proceeds accordingly
         *
         * - accepts idxs and row or col as args for sweep
         */
        void step_base_row(const T beta, const T q_max, int idx1, int idx2)
        {
            uniform_real_distribution<> get_q(-q_max, q_max);
            T old_mean = this->e_mean;
            T curr_H;

            this->steps++;

            gen_givens_row(idx1, idx2, get_q(gen));
            // rejection probability ~ e^(-beta * (H_f - H_i)), H_f > H_i
            // thresh can be > 1 if H_f < H_i, then always accept
            double reject_thresh = exp(beta * s_n * this->l1_delta);
            if (unit(gen) > reject_thresh)
            {
                this->rejects++;
            }
            else // below rejection threshold, accept
            {
                apply_givens_row(idx1, idx2);
            }
            // update statistics
            curr_H = this->H();
            this->e_mean = old_mean + (curr_H - old_mean) / steps;
            this->e_var_base += (curr_H - old_mean) * (curr_H - this->e_mean);
        }
        void step_base_col(const T beta, const T q_max, int idx1, int idx2)
        {
            uniform_real_distribution<> get_q(-q_max, q_max);
            T old_mean = this->e_mean;
            T curr_H;

            this->steps++;

            gen_givens_col(idx1, idx2, get_q(gen));
            // rejection probability ~ e^(-beta * (H_f - H_i)), H_f > H_i
            // thresh can be > 1 if H_f < H_i, then always accept
            double reject_thresh = exp(beta * s_n * this->l1_delta);
            if (unit(gen) > reject_thresh)
            {
                this->rejects++;
            }
            else // below rejection threshold, accept
            {
                apply_givens_col(idx1, idx2);
            }
            // update statistics
            curr_H = this->H();
            this->e_mean = old_mean + (curr_H - old_mean) / steps;
            this->e_var_base += (curr_H - old_mean) * (curr_H - this->e_mean);
        }
        void step(const T beta, const T q_max)
        {
            int idx1 = get_idx(gen);
            int idx2 = get_idx2(gen);
            if (idx2 >= idx1) idx2++; // ensure no duplicates
            if (row_or_col(gen) == 0)
            {
                step_base_row(beta, q_max, idx1, idx2);
            }
            else
            {
                step_base_col(beta, q_max, idx1, idx2);
            }
        }
        // for all
        void sweep(const T beta, const T q_max)
        {
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    step_base_row(beta, q_max, i, j);
                    step_base_col(beta, q_max, i, j);
                }
            }
        }

        // resets reject/step counter
        void reset_counts()
        {
            this->rejects = 0;
            this->steps = 0;
            this->e_mean = 0;
            this->e_var_base = 0;
        }
};
#endif
