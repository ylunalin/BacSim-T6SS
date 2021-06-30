#ifndef STRAIN_HH
#define STRAIN_HH
#include "common.hh"
#include "sim_params.hh"
#include "omp.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

struct strain {
	/** The bacterium type. */
    int type;
	/** The doubling time. */
    double t0;
	/** The average initial length. */
    double l0;
    /** A parameter that controls the division model. See Ariel Amir's paper. */
    double alpha;
    double sigma_t;
    double lag_time;
    double rate_on;

    strain() : type(-1), t0(-1), l0(-1), alpha(0), sigma_t(0) {}
    strain(int type_, sim_params *spars) :
        type(type_),
        t0(spars->double_times[type]),
        l0(spars->l0s[type]),
        alpha(spars->alphas[type]),
        sigma_t(t0 * spars->sigma_ts[type]),
        lag_time(spars->lag_times[type]),
        rate_on(spars->lag_times[type])
        {}

    double area0() const{
        return body_area_multiplier*l0 + cap_area;
    }

    double targ_length(const double lb) const{
        /* Linear expansion near l0 */
        //double talp = 2*alpha;
        //return talp * l0 + (2-talp) *lb + 2*r0;
        double tmp = 2* 2*pow((lb+r0), 1-alpha) * pow((l0+r0), alpha);
        return tmp;
    }

    double targ_length_w_noise(const double lb, double& td, const gsl_rng * rand) const{
        double  random_var = gauss_rnd(rand);
        td = t0* (1 + alpha* log((l0+r0)/(lb+r0)) / log(2))  + random_var;
        return targ_length(lb)*pow(2, random_var/t0);
    }

    double gauss_rnd(const gsl_rng * rand) const{
        // Here we draw a small gaussian
        // random variable, width is sigma_t
        return gsl_ran_gaussian(rand,sigma_t);
    }
};

#endif
