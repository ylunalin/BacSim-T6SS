/** A header file for the sim_params class, which contains
 * simulation parameters for the bac_sim class. Part of it
 * is adapted from fileinfo class in Chris's code. The
 * general structure of the code is similar to modularized dld.
 *
 * Author    : Y Luna Lin
 * Email     : y.lin2783@gmail.com
 * Date      : May 23, 2018
 */

#ifndef BAC_SIM_PARAMS
#define BAC_SIM_PARAMS

#include "common.hh"

/** Brief header file for a class that specifies simulation
 * parameters; will be particularly helpful when there are
 * many parameters we need to keep track and vary.
 */

class sim_params{
	public:
	/** All parameters and variables are dimensionless */
	/** Simulation domain config */
	bool x_prd, y_prd;
	int m, n, mn;
	double ax, bx, ay, by;

	/** Resource field config */
	bool diff_solve;
	int res_nx, res_ny;
	int diffusion_interval;
	double init_res_conc;
	double diffusion_const;


	/** Bacteria matrix config */
	// maximum radius to compute RDF for
	// default set to 20 body lenght (1 body length=2 b_rad)
	int bm_nx, bm_ny;
	int rdf_rmax;
    int chunk_size;

	/** Simulation duration config */
	/** Number of realizations */
	int nr;
    /** If to output bacteria data in all iterations */
    bool output_all;
	/** Number of frames */
	int frames;
	/** Threshold number of bacteria */
	int target_tb;
	/** Time step size */
	double dt;
	/** Total time duration */
	double T;
	/** Maximum allowed wall time for sims */
	double wall_time_max;

	/** Bacteria config */
	/** The type of initialization */
	int itype;
	/** Number of initial bacteria */
	int n_bac;
	/** Random number generator seed */
	int seed;
	/** The probability of type 1 */
	double ptype0;
	/** The size of the initilaization region. */
	double init_region_size;
	/** The secondary lengthscale in initial condition. */
	double init_sec_lengthscale;

    /** Initial length */
    double l0s[2];
    double lmax;
    /** Cell regulation factor alpha */
    double alphas[2];
	/** The growth rate of strains */
	double double_times[2];
	/** The lysis time of strains */
	double lysis_times[2];
    double min_db_t;
	/** The growth rate of strains, per hour */
	double base_growth_rates[2];
	/** The growth rate of strains, per hour */
	double growth_rates[2];
	/** The rates at which sheaths are produced, per hour */
    double sheath_rates[2];
    /** The average number of sheath of the initial Poisson dist */
    double init_mean_sh[2];
	/** Firing rate / division time unit, per hour*/
	double fire_rates[2];
	/** The size of noise in the interdivision time */
    double sigma_ts[2];
    double sigma_l;
    double sigma_a;

    /** Other internal model parameters*/
    double init_prob_on[2];
    double lag_times[2];
    double rate_on[2];

    /** The type of restriction on sheath production rate */
    restrict_type growth_factor;
    restrict_type sheath_factor;
	/** Growth yield */
	double gamma;
	/** A number dimensional number that measures
        ratio between elastic repulsive force and damping force. */
	double Edamp;
    double damp;
	/** Critial fraction of rep_force that would stop growth. */
	double stop_frac;
    /** Cost on growth rate due to sheath production. */
    double cost_coeff;

	// Experimental image configurations
	int width;
	int height;
	int dft_1D_profile_bin_size;
	bool negate_image;
	bool write_subplot;
	bool compute_dft;
	bool write_indv_dft;
	bool compute_rdf;
    bool compute_demixing;

	/** Output direcotry */
	char * dirname;


	/** Constructor that scans a config file. */
	sim_params(const char *config);

	/** Destructor */
	~sim_params(){
		delete [] dirname;
	}

	inline void transfer_basics(bool &xp, bool &yp, int &m_, int &n_, int &mn_, double &ax_, double &bx_, double &ay_, double &by_){
		xp = x_prd;
		yp = y_prd;
		m_ = m;
		n_ = n;
		mn_ = mn;
		ax_ = ax;
		bx_ = bx;
		ay_ = ay;
		by_ = by;
	}
	inline void change_rng_seed(int change){
		seed += change;
	}

	inline void print_report(){
		puts("#===============================================");
		printf("# Simulation parameters:\n"
               "# Rng seed: %d\n"
	       "# Periodicity: x: %d y: %d\n"
	       "# Nondimensional domain bounds: x: [%g, %g] y: [%g, %g]\n"
	       "# Num of blocks; x: %d y: %d\n"
	       "# Block sizes: sx: %g sy: %g\n#\n"
           "# Size of chunck in bacs: %d\n"

           "# Resource field number of nodes (%d %d)\n"
           "# nondimensional dx,dy = (%g %g)\n"
	       "# Solve diffusion or not %d\n"
           "# Nondimensional diffusion constant %g \n"
           "# Diffusion interval %d (nondimensional t = %g)\n"
           "# Nondimensional initial resource concentration %g\n#\n"

           "# Bacteria matrix field number of nodes (%d %d)\n"
           "# RDF maximum radius %d grid points\n#\n"

           "# Nondimensioanl bacterium setting\n"
           "#                    STRAIN 1  |  STRAIN 2\n"
	       "# Initial half length:   %g  %g\n"
	       "# Doubling time (T): %g  %g\n"
           "# Growth rate (T^-1): %g  %g\n"
           "# Noise strength in t_d : %g  %g\n"

           "# Sheath production rate (T^-1): %g  %g\n"
           "# Initial mean sheaths: %g  %g\n"
	       "# Firing rate (T^-1): %g  %g\n"
           "# Initial activated probability: %g %g\n"
           "# Lag time to begin (T): %g %g\n"
           "# Lysis times (T): %g %g\n"
           "# Activation rates: %g %g (T^-1)\n"
           "# Type of sheath production restriction: %s\n"
           "# Type of growth production restriction: %s\n"
	       "# Nondimensional growth yield ((g*s_h)/l_b): %g\n"
           "# Stopping force fraction %g\n"
	       "# Nondimensioanl number ratio of Elastic force to damping: %g\n"
           "# Nondimensional damping: %g\n"
           "# Unit repulsive force %g\n"
           "# Type 0 probability %g\n"
           "# Size of orientation noise %g\n"
           "# Size of length noise %g\n#\n"

	       "# Simulate to target number %d, then for T=%g, dt: %g, in %d frames\n"
           "# Repeat simulation %d times\n"
           "# Max wall time allowed to reach target %g s\n"
	       "# Initialize type: %s, number of bac: %d, region radius: %g\n#\n"

           "# Data analysis setting\n"
           "# Exp image dimension [%d x %d]\n"
           "# Negate image or not  %d\n"
           "# DFT 1D profile bin size %d\n"
           "# Write flags: subimages %d, individual 2D DFT %d\n"
           "# Compute flags: RDF %d, DFT %d, Demix %d\n#\n"

	       "# Output directory: %s\n",
		seed,
		x_prd, y_prd,
		ax, bx, ay, by,
		m, n,
		(bx-ax)/m, (by-ay)/n, chunk_size,

		res_nx, res_ny, x_prd?((bx-ax)/res_nx):((bx-ax)/(res_nx-1)),
		y_prd?((by-ay)/res_ny):((by-ay)/(res_ny-1)),
		diff_solve,
		diffusion_const, diffusion_interval,
		diffusion_interval*dt, init_res_conc,

		bm_nx, bm_ny, rdf_rmax,

		l0s[0], l0s[1],
		double_times[0], double_times[1],
        growth_rates[0], growth_rates[1],
        sigma_ts[0], sigma_ts[1],
        sheath_rates[0], sheath_rates[1],
        init_mean_sh[0], init_mean_sh[1],
		fire_rates[0], fire_rates[1],
        init_prob_on[0], init_prob_on[1],
        lag_times[0], lag_times[1],
        lysis_times[0], lysis_times[1],
        rate_on[0], rate_on[1],
        (sheath_factor==none)?"None":(sheath_factor==monod?"Monod only":(sheath_factor==crowd?"Crowding only":"Monod*Crowd")),
        (growth_factor==none)?"None":(growth_factor==monod?"Monod only":(growth_factor==crowd?"Crowding only":"Monod*Crowd")),
		gamma,
		stop_frac,
		Edamp,damp,Edamp * damp * lmax / min_db_t,ptype0,
        sigma_a, sigma_l,

		target_tb, T, dt, frames, nr, wall_time_max,
		itype==0?"Circle":(itype==1?"Square":(itype==2?"Line":"Unknown")), n_bac, init_region_size,

		width, height, negate_image, dft_1D_profile_bin_size, write_subplot, write_indv_dft,
		compute_rdf, compute_dft, compute_demixing,
		dirname);
		puts("#===============================================\n");
	}

	//private:
	/** Tests to see if two strings are equal.
	 * \param[in] p1 a pointer to the first string.
	 * \param[in] p2 a pointer to the second string.
	 * \return True if they are equal, false otherwise. */
	inline bool se(const char *p1,const char *p2) {
		return strcmp(p1,p2)==0;
	}
	/** Finds the next token in a string and interprets it as a
	 * double precision floating point number. If none is availble,
	 * it gives an error message.
	 * \param[in] ln the current line number. */
	inline double next_double(int ln) {
		return atof(next_token(ln));
	}
	/** Finds the next token in a string, interprets it as a double
	 * precision floating point number, and checks that there are
	 * no subsequent values.
	 * \param[in] ln the current line number. */
	inline double final_double(int ln) {
		double temp=next_double(ln);
		check_no_more(ln);
		return temp;
	}
	/** Finds the next token in a string and interprets it as an
	 * integer. If none is availble, it gives an error message.
	 * \param[in] ln the current line number. */
	inline int next_int(int ln) {
		return atoi(next_token(ln));
	}
	/** Finds the next token in a string, interprets it as an
	 * integer, and checks that there are no subsequent values.
	 * \param[in] ln the current line number. */
	inline int final_int(int ln) {
		int temp=next_int(ln);
		check_no_more(ln);
		return temp;
	}
	char* next_token(int ln);
	void invalid_error(const char *cs,int ln);
	void check_no_more(int ln);
};

#endif
