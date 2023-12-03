#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned int uint;

/* Function estimating the value of 'X' depending on the parameters. */
typedef void X_fct(double *X_mean, uint r, uint n, uint d, uint w, uint t,
                   uint S, double *rho_table);

/* Function modeling the sampling of a position. */
typedef void S_fct(double *pGarray, double *pBarray, uint n, uint d, uint t);

int max(int a, int b) { return a > b ? a : b; }

int min3(int a, int b, int c) {
    return a < b ? (a < c ? a : c) : (b < c ? b : c);
}

/* $lnbino(n, t) = \ln\binom{n}{t}$ */
double lnbino(uint n, uint t) {
    if ((t == 0) || (n == t))
        return 0.0;
    else
        return lgamma(n + 1) - lgamma(t + 1) - lgamma(n - t + 1);
}

double xlny(double x, double y) {
    if (x == 0.)
        return 0.;
    else
        return x * log(y);
}

/* Logarithm of the probability mass function of a binomial distribution:
 * $lnbinomial(n, k, p, q)) = \ln(\binom{n}{k} p^k q^{n-k})$ */
double lnbinomialpmf(int n, int k, double p, double q) {
    return lnbino(n, k) + xlny(k, p) + xlny(n - k, q);
}

/* $E_log(n, w, t, i) = \ln(\binom\binom{w}{i} {n - w}{t - i} / \binom{n}{t})$
 */
double E_log(uint n, uint w, uint t, uint i) {
    return lnbino(w, i) + lnbino(n - w, t - i) - lnbino(n, t);
}

/* X_avg = sum((l - 1) * E_l, l odd) */
double X_avg(uint r, uint n, uint w, uint t) {
    uint i;
    double x;

    /* E_log(n, w, t, i) decreases fast when 'i' varies.
     * For $i = 10$ it is very likely to be negligible. */
    for (x = 0, i = 1; (i <= w) && (i <= t); i += 2) {
        x += (i - 1) * exp(E_log(n, w, t, i));
    }

    return x * r;
}

/* Probability for a bit of the syndrome to be zero, knowing the syndrome
 * weight 'S' and 'X' */
double counters_C0(uint n, uint d, uint w, uint S, uint t, double x) {
    return ((w - 1) * S - x) / (n - t) / d;
}

/* Probability for a bit of the syndrome to be non-zero, knowing the syndrome
 * weight 'S' and 'X' */
double counters_C1(uint n, uint d, uint w, uint S, uint t, double x) {
    return (S + x) / t / d;
}

uint compute_threshold(uint r, uint n, uint d, uint w, uint S, uint t) {
    double p, q;

    double x = X_avg(r, n, w, t);
    p = counters_C0(n, d, w, S, t, x);
    q = counters_C1(n, d, w, S, t, x);

    uint threshold;
    if (q == 1.0 || p == 1.0 || p > q) {
        threshold = d;
    }
    else {
        threshold = d + 1;
        double diff = 0.;
        do {
            threshold--;
            diff = (-exp(lnbinomialpmf(d, threshold, p, 1. - p)) * (n - t) +
                    exp(lnbinomialpmf(d, threshold, q, 1. - q)) * t);
        } while (diff >= 0. && threshold > (d + 1) / 2);
        threshold = threshold < d ? (threshold + 1) : d;
    }

    return threshold;
}

/* The "locking" probability is the probability that no position among the 'n'
 * possible ones has a counter over the threshold.
 * The algorithm would then be "locked" and loop infinitely, leading to a
 * decoding failure. */
void lock(double *pGarray, double *pBarray, uint n, uint d, uint t, uint T,
          long double *p_lock) {
    long double pGlessT = 0.l;
    long double pBlessT = 0.l;
    for (uint x = 0; x < T; ++x) {
        pGlessT += pGarray[x];
        pBlessT += pBarray[x];
    }
    *p_lock = powl(pGlessT, n - t) * powl(pBlessT, t);
}

/* Fill *parray using a binomial probability mass function:
 * parray[i] = P[X = i]
 * where X follows a binomial distribution of parameters 'd' and 'p'. */
void binomial(double p, uint d, double *parray) {
    if (p >= 1) {
        for (uint x = 0; x < d; ++x) {
            parray[x] = 0.;
        }
        parray[d] = 1.;
    }
    else if (p <= 0) {
        for (uint x = 1; x <= d; ++x) {
            parray[x] = 0.;
        }
        parray[0] = 1.;
    }
    else {
        double q = 1. - p;

        for (uint x = 0; x <= d; ++x)
            parray[x] = exp(lnbinomialpmf(d, x, p, q));
    }
}

/* <-                n                ->
 * <-       w       ->
 * .-----------------------------------.
 * |1...............1|0...............0| <- row of the parity check matrix
 * '-----------------------------------'
 *
 * .-----------------------------------.
 * |1......1|0......0|1......1|0......0| <- error vector
 * '-----------------------------------'
 * <-  i  ->         <- t-i ->
 *
 * $rho_table[i] = \binom{w}{i} \binom{n - w}{t - i} / \binom{n}{t}$
 *
 * rho_table[i] is the probability that among all the positions in an equation
 * (a row of the parity check matrix) there are exactly 'i' errors.
 *
 * We assume that every row is chosen uniformly among the vectors of length 'n'
 * and Hamming weight 'w', and that the error vector is chosen uniformly among
 * the vectors of length 'n' and Hamming weight 't'.
 */
void rho(uint n, uint w, uint t, double *rho_table) {
    long double sum = 1.l;
    long double term = 1.l;
    rho_table[0] = 1.;
    for (uint i = 1; i <= w; ++i) {
        term *= (w - i + 1) * (t - i + 1);
        term /= i * (n - w - t + i);
        rho_table[i] = term;
        sum += term;
    }
    for (uint i = 0; i <= w; ++i) {
        rho_table[i] /= sum;
    }
}

long double sum_powers(long double p, uint K) {
    long double ret = 1.l;
    for (uint k = 0; k < K - 1; ++k)
        ret = 1.l + p * ret;
    return ret;
}

/* Assume that we chose 'T' as the threshold and knowing that the "good"
 * positions have counters distributed following *pGarray and the "bad"
 * positions' counters follow *pBarray.
 *
 * Modify *pGarray and *pBarray so that they contains the transition
 * probability i.e.
 * pBarray[i] = P[a position whose counter is 'i' is flipped and the position
 *               is an error];
 * pGarray[i] = P[a position whose counter is 'i' is flipped and the position
 *               is not an error].
 *
 * p_flip is the probability of flipping a position.
 * p_noflip is the probability of not flipping any position.
 */
void transitions(uint d, uint T, double *pGarray, double *pBarray,
                 long double *p_flip, long double *p_noflip) {
    *p_flip = 0.l;
    *p_noflip = 0.l;
    for (uint x = 0; x < T; ++x) {
        *p_noflip += pGarray[x] + pBarray[x];
        pGarray[x] = 0.;
        pBarray[x] = 0.;
    }
    for (uint x = T; x <= d; ++x) {
        *p_flip += pGarray[x] + pBarray[x];
    }
    *p_flip = *p_flip / (*p_flip + *p_noflip);
    *p_noflip = *p_noflip / (*p_flip + *p_noflip);
}

/* An unverified equation has an odd number of errors.
 *
 * Let us consider two quantities:
 * - the syndrome weight,
 * - and the sum of the counters of the errors.
 *
 * Say an equation has 'l' errors where 'l' is odd:
 * - it contributes to the syndrome by 1,
 * - it contributes to the sum of the counters of the errors by 'l'.
 * We call 'X' the difference between these two quantities.
 *
 * $X = (\sum_{l odd} (l - 1) * P[an equation has exactly l errors])
 *     /(\sum_{l odd} P[an equation has exactly l errors])$ */
void compute_X(uint n, uint w, uint t, double *X_mean, double *rho_table) {
    double term;
    long double den;

    den = 0.l;
    *X_mean = 0.;

    for (uint l = 1; l <= w; l += 2) {
        term = (l - 1) * rho_table[l];
        *X_mean += term;
        den += rho_table[l];
    }
    *X_mean = *X_mean / den;
}

/* Condition the distributions *pGarray and *pBarray to only the possible
 * counter values given 't' and 'S'. */
char condition_possible_counters(uint n, uint d, uint t, uint S,
                                 double *pGarray, double *pBarray) {
    long double denomG, denomB;

    denomG = 0.l;
    denomB = 0.l;
    /* For a "good" position (i.e. not an error) we have the following
     * inequations:
     * 0 <= sigma <= S
     * 0 <= S + d - 2 * sigma <= d * (t + 1)
     */
    uint iminG = max(0, ((int)S - (int)(d * t)) / 2);
    uint imaxG = min3(d, (S + d) / 2, S);
    /* For a "bad" position (i.e. an error) we have the following inequations:
     * 0 <= sigma <= S
     * 0 <= S + d - 2 * sigma <= d * (t - 1)
     */
    uint iminB = max(0, ((int)S - (int)d * ((int)t - 2)) / 2);
    uint imaxB = min3(d, (S + d) / 2, S);

    for (uint i = iminG; i <= imaxG; ++i) {
        denomG += pGarray[i];
    }
    for (uint i = iminB; i <= imaxB; ++i) {
        denomB += pBarray[i];
    }

    if (denomG == 0. || denomB == 0.) {
        return 0;
    }

    for (uint i = 0; i <= d; ++i) {
        if (i >= iminG && i <= imaxG)
            pGarray[i] /= denomG;
        else
            pGarray[i] = 0.;
    }
    for (uint i = 0; i <= d; ++i) {
        if (i >= iminB && i <= imaxB)
            pBarray[i] /= denomB;
        else
            pBarray[i] = 0.;
    }

    return 1;
}

/* Fill the array **dfr with decoding failure rates of an MDPC code:
 * dfr[t][S] = P[decoding a syndrome of weight 'S' corresponding to an error
 *              vector of weight 't' fails].
 *
 * The MDPC code considered has parameters 'r', 'n', 'd', 'w'.
 * Any state with t <= t_pass is considered to be sucessful for any syndrome
 * weight. Any state with t >= t_fail is considered a failure for any syndrome
 * weight.
 *
 * The decoding algorithm considered is the "step-by-step" decoder where the
 * position are sampled according to *sample_fct. The model for 'X' is given by
 * *x_fct.
 */
void fill_dfr_inf(uint r, uint n, uint d, uint w, uint t_pass, uint t_fail,
                  double (*dfr)[d * t_fail + 1], X_fct *x_fct,
                  S_fct *sample_fct) {
    long double(*p_lock)[d * t_fail + 1];
    long double(*p_flip)[d * t_fail + 1];
    long double(*p_noflip)[d * t_fail + 1];
    double *(*pGarray)[d * t_fail + 1];
    double *(*pBarray)[d * t_fail + 1];
    char(*reachable)[d * t_fail + 1];
    uint(*threshold)[d * t_fail + 1];

    p_lock = malloc((d * t_fail + 1) * t_fail * sizeof(long double));
    p_flip = malloc((d * t_fail + 1) * t_fail * sizeof(long double));
    p_noflip = malloc((d * t_fail + 1) * t_fail * sizeof(long double));
    pGarray = malloc((d * t_fail + 1) * t_fail * sizeof(double *));
    pBarray = malloc((d * t_fail + 1) * t_fail * sizeof(double *));
    reachable = calloc((d * t_fail + 1) * t_fail, sizeof(char));
    threshold = malloc((d * t_fail + 1) * t_fail * sizeof(uint));

#pragma omp parallel for
    for (uint S = 1; S <= r; ++S) {
        if (S > d * t_fail)
            continue;
        for (uint t = t_pass + 1; t < t_fail; ++t) {
            /* We know that the syndrome weight has the same parity as the
             * product d * t and that it cannot exceed d * t. We therefore
             * skip any other possibility. */
            if ((d * t % 2) != (S % 2) || S > d * t)
                continue;

            uint T2 = compute_threshold(r, n, d, w, S, t);
            threshold[t][S] = T2;

            double *rho_table = malloc((w + 1) * sizeof(double));
            double X;
            (*x_fct)(&X, r, n, d, w, t, S, rho_table);
            free(rho_table);

            double c0, c1;
            /* These arrays contain the probabilities for the counters, they
             * successively are used for:
             * - the distributions of the counters given the error weight 't'
             *   and the syndrome weight 'S';
             * - the distributions of the picked counters;
             * - the distributions of the flipped counters.
             *
             * B is for "bad" positions (i.e. errors)
             * G is for "good" positions */
            pBarray[t][S] = malloc((d + 1) * sizeof(double));
            pGarray[t][S] = malloc((d + 1) * sizeof(double));

            /* First, compute the counter distributions:
             * pBarray[i] = P[sigma = i knowing that the position is an error]
             * pGarray[i] = P[sigma = i knowing that the position is not an
             * error] where sigma is the counter of a position (the number of
             * unverified equations it is involved in).
             *
             * We assume they follow a binomial distributions of parameters:
             * - d and c0 for "good" positions;
             * - d and c1 for "bad" positions.
             * */
            c0 = ((w - 1) * S - X) / (d * (n - t));
            c1 = (S + X) / (d * t);
            binomial(c1, d, pBarray[t][S]);
            binomial(c0, d, pGarray[t][S]);

            dfr[t][S] = 1.;
            /* Some counters are not possible given 't' and 'S' (for example
             * when only one error is left, its only counter value possible is
             * 'd').
             *
             * This function conditions the probability knowing the only
             * possible values possible). If no counter value is possible, we
             * skip this state. */
            if (!condition_possible_counters(n, d, t, S, pGarray[t][S],
                                             pBarray[t][S])) {
                free(pBarray[t][S]);
                free(pGarray[t][S]);
                continue;
            }
            else {
                reachable[t][S] = 1;
            }

            /* Compute the "locking" probability (i.e. the probability that no
             * position has a counter over the threshold). */
            lock(pGarray[t][S], pBarray[t][S], n, d, t, threshold[t][S],
                 &p_lock[t][S]);
            /* Compute the counters probabilities after the sampling method. */
            (*sample_fct)(pGarray[t][S], pBarray[t][S], n, d, t);
            /* Finally consider the threshold.
             * Compute the probability of flipping a position depending on the
             * fact that it is an error or not and on its counter value. The
             * total probability of flipping a position is 'p_flip'. The
             * probability to flip no position is 'p_noflip'. */
            transitions(d, threshold[t][S], pGarray[t][S], pBarray[t][S],
                        &p_flip[t][S], &p_noflip[t][S]);
        }
    }

    /* Consider any state with a zero syndrome weight and a nonzero error
     * weight to be a failure. */
    for (uint t = t_pass + 1; t < t_fail; ++t) {
        if ((d * t % 2) != 0)
            continue;
        dfr[t][0] = 1.;
    }

    for (uint S = 1; S <= d * t_fail; ++S) {
#pragma omp parallel for
        for (uint t = max(t_pass + 1, S / d); t < t_fail; ++t) {
            if (!reachable[t][S])
                continue;

            double *pGarray_current = pGarray[t][S];
            double *pBarray_current = pBarray[t][S];
            long double p_lock_current = p_lock[t][S];
            long double p_flip_current = p_flip[t][S];
            long double p_noflip_current = p_noflip[t][S];
            long double dfr_current;
            int flip = 0;
            /* Probabilities over 1 can arise given rounding errors and
             * other imprecisions, we change them into 1 to avoid
             * propagation of these computation errors. */
            if (p_noflip_current < 1.) {
                dfr_current = 0.l;
                uint T = threshold[t][S];
                for (int Sigma = T; Sigma <= d; ++Sigma) {
                    if (pGarray_current[Sigma] != 0.) {
                        flip |= 1;
                        if (t + 1 < t_fail)
                            dfr_current += pGarray_current[Sigma] *
                                           dfr[t + 1][S + d - 2 * Sigma];
                        else {
                            /* Assume this transition leads to a failure. */
                            if (pGarray_current[Sigma] != 0.)
                                dfr_current += pGarray_current[Sigma];
                        }
                    }
                    if (pBarray_current[Sigma] != 0.) {
                        flip |= 1;
                        if (t - 1 > t_pass) {
                            dfr_current += pBarray_current[Sigma] *
                                           dfr[t - 1][S + d - 2 * Sigma];
                        }
                    }
                }
            }
            if (!flip)
                dfr[t][S] = 1.;
            else {
                /* We suppose that the algorithm can perform an infinite number
                 * of iterations until it finds a postion to flip. */
                dfr[t][S] = p_lock_current + (1.l - p_lock_current) *
                                                 (dfr_current / p_flip_current);
                dfr[t][S] = (dfr[t][S] > 1.) ? 1. : dfr[t][S];
            }
        }
    }
    for (uint S = 1; S <= d * t_fail; ++S) {
#pragma omp parallel for
        for (uint t = t_pass + 1; t < t_fail; ++t) {
            if (reachable[t][S]) {
                free(pGarray[t][S]);
                free(pBarray[t][S]);
            }
        }
    }
    free(p_lock);
    free(p_flip);
    free(p_noflip);
    free(pGarray);
    free(pBarray);
    free(reachable);
    free(threshold);
}

/* Assume 'X' is always 0.
 * This is a quite pessimistic hypothesis. */
void noX(double *X, uint r, uint n, uint d, uint w, uint t, uint S,
         double *rho_table) {
    *X = 0.;
}

/* Assume 'X' depends only on the Hamming weight of the error vector 't'. */
void Xt(double *X, uint r, uint n, uint d, uint w, uint t, uint S,
        double *rho_table) {
    rho(n, w, t, rho_table);
    compute_X(n, w, t, X, rho_table);
    *X *= r;
}

/* Assume 'X' depends on the Hamming weight of the error vector 't' and the
 * weight of the syndrome 'S'. */
void XSt(double *X, uint r, uint n, uint d, uint w, uint t, uint S,
         double *rho_table) {
    rho(n, w, t, rho_table);
    compute_X(n, w, t, X, rho_table);
    *X *= S;
}

/* Assume positions are sampled uniformly.
 *
 * Modify *pGarray and *pBarray so that they follow:
 * pBarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               an error];
 * pGarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               not an error]. */
void uniform(double *pGarray, double *pBarray, uint n, uint d, uint t) {
    for (uint x = 0; x <= d; ++x) {
        pGarray[x] *= (double)(n - t) / n;
        pBarray[x] *= (double)t / n;
    }
}

/* Assume positions are chosen by first picking a random unverified equation
 * then picking a random postion in this equation.
 *
 * Modify *pGarray and *pBarray so that they follow:
 * pBarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               an error];
 * pGarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               not an error]. */
void single(double *pGarray, double *pBarray, uint n, uint d, uint t) {
    long double total = 0.l;
    for (uint x = 0; x <= d; ++x) {
        pGarray[x] *= x * (n - t);
        pBarray[x] *= x * t;
        total += pGarray[x] + pBarray[x];
    }
    if (total > 0.) {
        for (uint x = 0; x <= d; ++x) {
            pGarray[x] /= total;
            pBarray[x] /= total;
        }
    }
}

/* Assume positions are chosen by first picking two random unverified equations
 * then picking a random postion in this two equations.
 *
 * Modify *pGarray and *pBarray so that they follow:
 * pBarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               an error];
 * pGarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               not an error]. */
void pair(double *pGarray, double *pBarray, uint n, uint d, uint t) {
    long double total = 0.l;
    for (uint x = 0; x <= d; ++x) {
        pGarray[x] *= x * (x - 1) * (n - t);
        pBarray[x] *= x * (x - 1) * t;
        total += pGarray[x] + pBarray[x];
    }
    if (total > 0.) {
        for (uint x = 0; x <= d; ++x) {
            pGarray[x] /= total;
            pBarray[x] /= total;
        }
    }
}

/* Assume positions are chosen by first picking three random unverified
 * equations then picking a random postion in this three equations.
 *
 * Modify *pGarray and *pBarray so that they follow:
 * pBarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               an error];
 * pGarray[i] = P[a position whose counter is 'i' is picked and the position is
 *               not an error]. */
void triple(double *pGarray, double *pBarray, uint n, uint d, uint t) {
    long double total = 0.l;
    for (uint x = 0; x <= d; ++x) {
        pGarray[x] *= x * (x - 1) * (x - 2) * (n - t);
        pBarray[x] *= x * (x - 1) * (x - 2) * t;
        total += pGarray[x] + pBarray[x];
    }
    if (total > 0.) {
        for (uint x = 0; x <= d; ++x) {
            pGarray[x] /= total;
            pBarray[x] /= total;
        }
    }
}

void print_help(char *arg0) {
    printf(
        "Usage: %s r d t_pass t_fail [noX|XSt|Xt] [uniform|single|pair|triple]\n"
        "\n"
        "Arguments:\n"
        "  r                     - Block length of the QC-MDPC code\n"
        "  d                     - Block weight of the QC-MDPC code\n"
        "  t_pass                - Threshold value; decoder is assumed to always succeed below this\n"
        "  t_fail                - Threshold value; decoder is assumed to always fail above this\n"
        "  [noX|XSt|Xt]          - Counters modelling choices:\n"
        "                           noX: assume X is always 0\n"
        "                           XSt: approximate X depending on S and t\n"
        "                           Xt: only consider t\n"
        "  [uniform|single|pair|triple]\n"
        "                        - Position choosing strategies:\n"
        "                           uniform: choose a position uniformly at random\n"
        "                           single: choose a random position involved in a random unverified equation\n"
        "                           pair: choose a random position involved in two random unverified equations\n"
        "                           triple: choose a random position involved in three random unverified equations\n"
        "\n"
        "Example: %s 12323 71 5 150 XSt single\n",
        arg0, arg0);
    exit(1);
}

int main(int argc, char *argv[]) {
    if (argc <= 6) {
        print_help(argv[0]);
    }
    uint r = atoi(argv[1]);
    uint n = 2 * r;
    uint d = atoi(argv[2]);
    uint w = 2 * d;
    uint t_pass = atoi(argv[3]);
    uint t_fail = atoi(argv[4]);
    X_fct *x_fct = NULL;
    S_fct *sample_fct = NULL;

    if (!strcmp(argv[5], "noX")) {
        x_fct = noX;
    }
    else if (!strcmp(argv[5], "XSt")) {
        x_fct = XSt;
    }
    else if (!strcmp(argv[5], "Xt")) {
        x_fct = Xt;
    }
    else {
        print_help(argv[0]);
    }

    if (!strcmp(argv[6], "uniform")) {
        sample_fct = uniform;
    }
    else if (!strcmp(argv[6], "single")) {
        sample_fct = single;
    }
    else if (!strcmp(argv[6], "pair")) {
        sample_fct = pair;
    }
    else if (!strcmp(argv[6], "triple")) {
        sample_fct = triple;
    }
    else {
        print_help(argv[0]);
    }

    double(*dfr)[d * t_fail + 1] =
        malloc((d * t_fail + 1) * t_fail * sizeof(double));

    fill_dfr_inf(r, n, d, w, t_pass, t_fail, dfr, x_fct, sample_fct);

    for (uint t = t_pass + 1; t < t_fail; ++t) {
        for (uint S = (d * t) % 2; S <= d * t && S <= r; S += 2) {
            printf("%d %d %.50g\n", t, S, dfr[t][S]);
        }
    }

    free(dfr);

    exit(EXIT_SUCCESS);
}
