/** Implementation of sim_params class.
 *
 * Author    : Y Luna Lin
 * Email     : y.lin2783@gmail.com
 * Date      : May 23, 2018
 */

#include "common.hh"
#include "sim_params.hh"

sim_params::sim_params(const char *fn):
	// Populate the parameters with default values
	// In case the config file does not specify
	x_prd(1), y_prd(1), m(25), n(25), mn(m*n),
	ax(-50.), bx(50.), ay(-50), by(50.),
	diff_solve(true),
	res_nx(100), res_ny(100), init_res_conc(1.),
	diffusion_const(1),
	bm_nx((int) ((bx-ax) / (DEFAULT_SIM_PARS::d_b_rad*1.725) + 0.5)),
	bm_ny((int) ((by-ay) / (DEFAULT_SIM_PARS::d_b_rad*1.725) + 0.5)),
	rdf_rmax((int) (bm_nx/2)), chunk_size(400),
	nr(1), output_all(0), frames(300), target_tb(4000),
	dt(1.), T(50.), wall_time_max(18000.),
	itype(1), n_bac(25), seed(1234), ptype0(0.5),
	init_region_size(50.), init_sec_lengthscale(0.5*init_region_size),

	l0s{DEFAULT_SIM_PARS::d_init_l, DEFAULT_SIM_PARS::d_init_l},
    alphas{1., 1.},
    double_times{0.75, 0.75},
    lysis_times{0.75, 0.75},
    base_growth_rates{0.924, 0.924},
    growth_rates{0.924, 0.924},
    sheath_rates{1., 1.},
    init_mean_sh{1., 1.},
	fire_rates{1.,1.},
    sigma_ts{0.02, 0.02},
    sigma_l(0.05),
    sigma_a(0.01),
    init_prob_on{0.,0.},
    lag_times{0,0},
    rate_on{10,10},

    growth_factor(monod),
    sheath_factor(monod),
    gamma(DEFAULT_SIM_PARS::d_gamma),
	Edamp(DEFAULT_SIM_PARS::d_Edamp),
    damp(DEFAULT_SIM_PARS::d_damp),
	stop_frac(DEFAULT_SIM_PARS::d_stop_frac),
    cost_coeff(DEFAULT_SIM_PARS::cost_coeff),

	width(1024), height(1344), dft_1D_profile_bin_size(1),
	negate_image(false), write_subplot(false),
	compute_dft(false), write_indv_dft(false),
    compute_rdf(false), compute_demixing(false)
	{
	int i, ln=1, len = strlen(fn);
	char buf[DEFAULT_SIM_PARS::buf_size], *bp;
	if(len<4 || fn[len-4]!='.' || fn[len-3]!='c' || fn[len-2]!='f' || fn[len-1]!='g'){
		fatal_error("sim_params::sim_params(const char*): Filename must end in '.cfg'\n", 1);
	}

	// Use the root of config filename as output directory name
	dirname = new char[len+1];
	for(i=0;i<len-3;i++) dirname[i]=fn[i];
	dirname[len-3]='o';
	dirname[len-2]='u';
	dirname[len-1]='t';
    dirname[len]='\0';
	FILE *f=safe_fopen(fn,"r");
    if(f==NULL) {
        exit(1);
    }
	while(!feof(f)){
		// Reads in config file one line at a time
		if(fgets(buf, DEFAULT_SIM_PARS::buf_size, f)==NULL) break;

		// Replace comment character with a null character, skip commented lines
		bp=buf;
		while((*bp)!=0){
			if(*bp=='#') {*bp=0; break;}
			bp++;
		}

		// Search for a keyword, skip if none found
		bp = strtok(buf, " \t\n");
		if(bp==NULL){
			ln++;	continue;
		}
		// Compare keyword with known keywords
		if(se(bp, "x_periodic")){
			x_prd = final_int(ln);
		} else if(se(bp, "y_periodic")){
			y_prd = final_int(ln);
		} else if(se(bp, "xy_bounds")){
			ax = next_double(ln);
			bx = next_double(ln);
			ay = next_double(ln);
			by = final_double(ln);
		} else if(se(bp, "grid_dims")){
			m = next_int(ln);
			n = final_int(ln);
			mn = m*n;
        } else if(se(bp, "init_lens")){
            l0s[0] = next_double(ln);
            l0s[1] = final_double(ln);
        } else if(se(bp, "alphas")){
            alphas[0] = next_double(ln);
            alphas[1] = final_double(ln);
		} else if(se(bp, "double_times")){
			double_times[0] = next_double(ln);
			double_times[1] = final_double(ln);
		} else if(se(bp, "lysis_times")){
			lysis_times[0] = next_double(ln);
			lysis_times[1] = final_double(ln);
		} else if(se(bp, "sheath_rates")){
			sheath_rates[0] = next_double(ln);
			sheath_rates[1] = final_double(ln);
		} else if(se(bp, "init_mean_sheaths")){
			init_mean_sh[0] = next_double(ln);
			init_mean_sh[1] = final_double(ln);
		} else if(se(bp, "fire_rates")){
			fire_rates[0] = next_double(ln);
			fire_rates[1] = final_double(ln);
        } else if(se(bp, "sigma_l")){
            sigma_l = next_double(ln);
        } else if(se(bp, "sigma_a")){
            sigma_a = next_double(ln);
		} else if(se(bp, "sigma_ts")){
			sigma_ts[0]= next_double(ln);
			sigma_ts[1]= final_double(ln);
		} else if(se(bp, "init_prob_on")){
			init_prob_on[0]= next_double(ln);
			init_prob_on[1]= final_double(ln);
		} else if(se(bp, "lag_times")){
			lag_times[0]= next_double(ln);
			lag_times[1]= final_double(ln);
		} else if(se(bp, "rate_on")){
			rate_on[0]= next_double(ln);
			rate_on[1]= final_double(ln);
        } else if(se(bp, "growth_factor")){
            growth_factor = (restrict_type)final_int(ln);
        } else if(se(bp, "sheath_factor")){
            sheath_factor = (restrict_type)final_int(ln);
		} else if(se(bp, "stop_frac")){
			stop_frac = final_double(ln);
		} else if(se(bp, "width")){
			width= final_int(ln);
		} else if(se(bp, "height")){
			height= final_int(ln);
		} else if(se(bp, "dft_1D_profile_bin_size")){
			dft_1D_profile_bin_size = final_int(ln);
		} else if(se(bp, "negate_image")){
			negate_image= final_int(ln);
		} else if(se(bp, "compute_dft")){
			compute_dft= final_int(ln);
		} else if(se(bp, "write_indv_dft")){
			write_indv_dft = final_int(ln);
        } else if(se(bp, "compute_demixing")){
            compute_demixing= final_int(ln);
		} else if(se(bp, "compute_rdf")){
			compute_rdf= final_int(ln);
		} else if(se(bp, "write_subplot")){
			write_subplot = final_int(ln);
		} else if(se(bp, "gamma")){
			gamma = final_double(ln);
		} else if(se(bp, "init_type")){
			int it= final_int(ln);
			if(it>=0 && it<=4) itype =it;
			else{
				printf("Initialization type %d not recognized, set to default %d\n", it, itype);
			}
		} else if(se(bp, "Edamp")){
			Edamp = final_double(ln);
		} else if(se(bp, "damp")){
			damp = final_double(ln);
		} else if(se(bp, "init_num_bac")){
			n_bac = final_int(ln);
		} else if(se(bp, "seed")){
			seed = final_int(ln);
		} else if(se(bp, "prob_type0")){
			ptype0 = final_double(ln);
		} else if(se(bp, "init_region_size")){
			init_region_size = final_double(ln);
		} else if(se(bp, "init_sec_lengthscale")){
			init_sec_lengthscale = final_double(ln);
		} else if(se(bp, "dt")){
			dt = final_double(ln);
		} else if(se(bp, "T")){
			T = final_double(ln);
        } else if(se(bp, "output_all")){
            output_all = final_int(ln);
		} else if(se(bp, "realization")){
			nr = final_int(ln);
		} else if(se(bp, "frames")){
			frames = final_int(ln);
		} else if(se(bp, "target_tb")){
			target_tb = final_int(ln);
		} else if(se(bp, "wall_time_max")){
			wall_time_max = final_double(ln);
		} else if(se(bp, "res_nx")){
			res_nx = final_int(ln);
		} else if(se(bp, "res_ny")){
			res_ny = final_int(ln);
		} else if(se(bp, "init_res_conc")){
			init_res_conc = final_double(ln);
		} else if(se(bp, "diff_solve")){
			diff_solve = (final_int(ln)>0);
		} else if(se(bp, "diffusion_const")){
			diffusion_const = final_double(ln);
		} else if(se(bp, "cost_coeff")){
			cost_coeff = final_double(ln);
		} else if(se(bp, "rdf_rmax")){
			rdf_rmax = final_int(ln);
		} else if(se(bp, "chunk_size")){
			chunk_size = final_int(ln);
		} else {
			printf("sim_params::sim_params(const char*): Unrecognized keyword '%c' at line %d of file %s\n", *bp, ln, fn);
			exit(1);
		}
		ln++;

	}

    // !! Some parameters are dependent on other
    // So we recalculate those one, in case the dependencies have changed.

	// Check if the initialization region fits within the simulation domain
	double xmin = fabs(ax)<fabs(bx)?fabs(ax):fabs(bx);
	double ymin = fabs(ay)<fabs(by)?fabs(ay):fabs(by);
	xmin = xmin<ymin?xmin:ymin;
	if(init_region_size>xmin) init_region_size = xmin;
	if(init_sec_lengthscale>xmin) init_sec_lengthscale = xmin;

	bm_nx=(int) ((bx-ax) / (DEFAULT_SIM_PARS::d_b_rad*1.725) + 0.5);
	bm_ny=(int) ((by-ay) / (DEFAULT_SIM_PARS::d_b_rad*1.725) + 0.5);
	rdf_rmax=(int) (bm_nx/2);
    if(chunk_size<200) chunk_size=200;

    lmax = l0s[0]>l0s[1]?l0s[0]:l0s[1];
    min_db_t = double_times[0]<double_times[1]?double_times[0]:double_times[1];
    for(int i=0;i<2;i++) {
            base_growth_rates [i] = log(2) / double_times[i];
            growth_rates[i] = base_growth_rates[i] - cost_coeff * sheath_rates[i];
    }
    double max_growth_rate =  (base_growth_rates[0]>base_growth_rates[1])?base_growth_rates[0]:base_growth_rates[1];
    double min_damp = max_growth_rate*(2*lmax + 1.33333333*r0) * cap_area;
    if(damp<min_damp) damp = min_damp;

    double lmin = l0s[0]>l0s[1]?l0s[1]:l0s[0];
    double mass_min = (2*lmin + 1.33333333*r0) * cap_area;
    double ideal_dt = mass_min / (Edamp * damp);
	if(dt>ideal_dt) dt = ideal_dt;

	// Compute a rough estimate of number of timesteps to take per each diffusion solve
	double dx = x_prd?((bx-ax)/res_nx):((bx-ax)/(res_nx-1));
	double dy = y_prd?((by-ay)/res_ny):((by-ay)/(res_ny-1));
	double dh = (dx>dy)?dy:dx;
	double tmp = (dh*dh/diffusion_const) / dt + 0.5;
	if(tmp>(20./dt)) tmp=(20./dt);
	diffusion_interval =  int(tmp);

	print_report();
}

/** Finds the next token in a string, and if none is availble, gives an error
 * message.
 * \param[in] ln the current line number. */
char* sim_params::next_token(int ln) {
    char *temp=strtok(NULL," \t\n");
    if(temp==NULL) {
        fprintf(stderr,"sim_params::Not enough arguments at input line %d\n",ln);
        exit(1);
    }
    return temp;
}

/** Prints a message about an invalid quantity and exits.
 * \param[in] cs the invalid quantity.
 * \param[in] ln the current line number. */
void sim_params::invalid_error(const char* cs,int ln) {
    fprintf(stderr,"sim_params::Invalid %s at line %d",cs,ln);
    exit(1);
}

/** Checks that there are no subsequent values.
 * \param[in] ln the current line number. */
void sim_params::check_no_more(int ln) {
    if(strtok(NULL," \t\n")!=NULL) {
        fprintf(stderr,"sim_params::Too many arguments at input line %d\n",ln);
        exit(1);
    }
}
