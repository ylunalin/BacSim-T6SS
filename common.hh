#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

const double pi=3.1415926535897932384626433832795;
const double twopi=2*3.1415926535897932384626433832795;
const double halfpi=0.5*3.1415926535897932384626433832795;

namespace DEFAULT_SIM_PARS{

    const int buf_size=1024;

    const double d_init_l = 1;

    /** Default growth rate (equal) for two strains.
        ~0.75 T doubling time. */
    const double d_gr = 0.925;

    /** Default kill rate (equal) for two strains. */
    const double d_kr = 0.1;

    /** Default growth yield. */
    const double d_gamma = 2.5;

    /** Default stopping_force percentage. */
    const double d_stop_frac = 0.1;

    /** Default radius of cells */
    const double d_b_rad=1;

    /** The nondimensional number measuring ratio
     * between Elastic repulsive force and damping force. */
    const double d_Edamp=10;

    const double d_damp=d_gr*(1.3333* d_b_rad + 2*d_init_l)* (pi*d_b_rad*d_b_rad);

    const double cost_coeff=0.001;

}

enum restrict_type: int{
    none = 0,
    monod = 1,
    crowd = 2,
    monod_n_crd = 3
};

const double r0 = DEFAULT_SIM_PARS::d_b_rad;
const double cap_area = r0*r0*pi;
const double body_area_multiplier = 4*r0;

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
double argument(double x,double y);

#endif
