#include "bac_sim.hh"

/** Constructs a bac_sim class, allocating a block structure to store the
 * bacteria positions and opening a diagnostic output file.
 * \param[in] sp the pointer to simulation parameters object
 * \param[in] (ax_,bx_) the x dimensions of the simulation region.
 * \param[in] (ay_,by_) the y dimensions of the simulation region.
 * \param[in] (x_prd_,y_prd_) the periodicity.
 * \param[in] (m,n) the number of blocks to divide the simulation region into.
 * \param[in] filename the directory name in which to store the output. */
bac_sim::bac_sim(sim_params *sp):
	spars(sp),
    x_prd(spars->x_prd), y_prd(spars->y_prd),
    m(spars->m), n(spars->n), mn(m*n),
	ax(spars->ax), bx(spars->bx), ay(spars->ay), by(spars->by),
	dx((bx-ax)/m), dy((by-ay)/n), xsp(1./dx),
	ysp(1./dy), filename(sp->dirname),
	n_bact(0), f_count(0), nr(spars->nr),
    time(0),
	co(new int[mn]), gh(new int[mn]), mem(new int[mn]),
	ba(new bact*[mn]),
    bac_types{strain(0,spars), strain(1,spars)},
    compt_times{0,0,0},
	// Private members:
	internal_sp(false), jmin(n), jmax(0), threads(1), iter(0), tot_l(0),
#ifdef _OPENMP
	//max_threads(1),
	max_threads(omp_get_max_threads()>20?20:omp_get_max_threads()),
#else
	max_threads(1),
#endif
    two_b_rad(DEFAULT_SIM_PARS::d_b_rad*2), two_b_rad_rsq(two_b_rad*two_b_rad),
	l_thresh(two_b_rad+2*(spars->l0s[0] + spars->l0s[1])),
	F_r(spars->Edamp * spars->damp * spars->lmax / spars->min_db_t),
	gamma(spars->gamma), gsp(1./gamma),
	stop_force(spars->stop_frac*F_r),
	buf(new char[256]),
	randos(new gsl_rng*[max_threads]) {
#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#endif
	int i;
	// Seed the random number generator to a specific value
	srand(sp->seed);

	// Initialize the simulation blocks
	for(i=0;i<mn;i++) {
		co[i]=0;
		mem[i]=init_region_memory;
		ba[i]=new bact[init_region_memory];
	}

	// Initialize the thread-base random number generators
	// One per thread
	for(i=0;i<max_threads;i++) {
		randos[i]=gsl_rng_alloc(gsl_rng_taus);
	}

	set_rnd(sp->seed);
    printf("\n# Maximum number of threads %d\n", max_threads);
    printf("# Elastic k = %g, sqrt(k/m)^-1 = %g, dt = %g\n", F_r, 1./sqrt(F_r/10), spars->dt);
    printf("# Total area available for growth %g\n", (bx-ax)*(by-ay));
    printf("# Stopping force %g\n\n", stop_force);

}

/** The class destructor closes the diagnostic file and frees the dynamically
 * allocated memory. */
bac_sim::~bac_sim() {
	for(int i=0;i<max_threads;i++) gsl_rng_free(randos[i]);
	delete [] randos;
	delete [] buf;
	if(internal_sp) delete spars;
	for(int i=mn-1;i>=0;i--) delete [] ba[i];
	delete [] ba;
	// delete [] w;
	delete [] mem;
	delete [] gh;
	delete [] co;
}

void bac_sim::reset(int iter_, int seed_){
	iter = iter_;
	set_rnd(seed_);
	tot_l = 0;

    if(n_bact!=0) {
        for(int i=0;i<mn;i++) {
            bact * nba = new bact[init_region_memory];
            delete []  ba[i];
            ba[i] = nba;
            co[i] = 0; gh[i] = 0; mem[i] = init_region_memory;
        }
    }
	n_bact = 0; f_count = 0; time = 0;
}

/** Call the solve(...) function using info in sim_params object */
void bac_sim::solve(bool output, bool verbose, double *data){
    if(internal_sp){
        printf("!!!!!!! Warning: There is an internal sim_params object with default options."
               " Are you sure you want to use default duration, frames params?\n");
    }
    double addition_T = spars->T - time;
    int addition_frames = spars->frames - f_count;

    if(addition_T < 0) addition_T = 0.;
    if(addition_frames < 1) addition_frames=1;
    solve(addition_T, addition_frames, spars->dt, output, verbose, data);

}

/** Call the solve(...) function using info in sim_params object until a target is reached. */
void bac_sim::solve(bool output, bool verbose){
    if(internal_sp){
        printf("!!!!!!! Warning: There is an internal sim_params object with default options."
               " Are you sure you want to use default duration, frames params?\n");
    }
    solve(spars->dt, spars->wall_time_max, spars->target_tb, output, verbose);
}


/** Carries out a bacterial growth simulation
 * \param[in] duration the time interval over which to simulate.
 * \param[in] frames the number of frames to store. */
void bac_sim::solve(double duration,int frames,double dt_max,bool output, bool verbose, double *data) {
    if(duration <= 0. || frames <1) return;

	double t_start=time,time_interval=duration/frames,target_time,t,t0,t2;
	int l=0,tb=0,num_st1=0,num_st2=0;
	double ap1=0, ap2=0;
    double sh1=0, sh2=0;

	// Output the initial fields and record initial time
	//tb = total_bacteria();
    double total_area=0;
    // number of cells with crd factor = 0
    int nzcf = 0;
	specific_bac_count(num_st1, num_st2, nzcf, ap1, ap2, total_area, sh1, sh2);
	tb = num_st1+num_st2;
	if(verbose){
		printf("# Frame\t time\t nbact (zero growth)\t strain0(red)\t strain1(blue)\t active0\t active1\t sheath1\t sheath2\t occupied area (perc of tot) \t threads\t calc_force time\t int_and_div time\t remap time\n");
		printf("# %d\t%.04g\t%d (%d) \t%d\t%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.06g\t%.04g\t%d\t --\t--\t--\n", f_count, time, tb, nzcf, num_st1, num_st2, ap1, ap2, sh1, sh2, total_area, total_area*100/((bx-ax)*(by-ay)), threads);
	}
    int nfields = 14;
	// storing data to external structure:
	data[0] = time;
	data[1] += double(tb)/double(nr);
	data[2] += double(num_st1) / double(nr);
	data[3] += double(num_st2)/double(nr);
	data[4] += double(num_st1)*double(num_st1)/double(nr);
	data[5] += double(num_st2)*double(num_st2)/double(nr);
	data[6] += double(ap1) / double(nr);
	data[7] += double(ap2)/double(nr);
	data[8] += double(ap1)*double(ap1)/double(nr);
	data[9] += double(ap2)*double(ap2)/double(nr);
	data[10] += double(sh1) / double(nr);
	data[11] += double(sh2)/double(nr);
	data[12] += double(sh1)*double(sh1)/double(nr);
	data[13] += double(sh2)*double(sh2)/double(nr);

	if(f_count==0 && output) {
		write(0);
	}

	t0=wtime();t=t0;
	for(int k=1;k<=frames;k++) {

		// Compute the target time to the next output frame
		target_time=t_start+time_interval*k;l=0;

		// Calculate threads
		threads=n_bact/spars->chunk_size+1;
		if(threads>max_threads) threads=max_threads;

		// Carry out simulation step using the regular timestep until
		// within range of the target time
		while(time+dt_max*(1+1e-8)<target_time) {
			l++;
			step_forward(dt_max);
		}

		// Carry out a final simulation step, using exactly the right
		// time step to reach the target time
		step_forward(target_time-time);
		l++;
		// Print diagnostic information
		t2=wtime();
		specific_bac_count(num_st1, num_st2, nzcf, ap1, ap2, total_area, sh1, sh2);
		tb = num_st1+num_st2;
		if(verbose){
            printf("# %d\t%.04g\t%d (%d) \t%d\t%d \t %.04g \t %.04g \t%.04g \t%.04g \t%.06g \t%.02g\t %d \t%.04g \t%.04g \t%.04g\n", k+f_count, time, tb, nzcf, num_st1, num_st2, ap1, ap2, sh1, sh2, total_area, total_area*100/((bx-ax)*(by-ay)), threads, compt_times[0], compt_times[1], compt_times[2]);
	}
		// storing data to external structure:
		data[k*nfields] = time;
		data[k*nfields+1] += double(tb)/double(nr);
		data[k*nfields+2] += double(num_st1) / double(nr);
		data[k*nfields+3] += double(num_st2)/double(nr);
		data[k*nfields+4] += double(num_st1)*double(num_st1)/double(nr);
		data[k*nfields+5] += double(num_st2)*double(num_st2)/double(nr);
		data[k*nfields+6] += double(ap1) / double(nr);
		data[k*nfields+7] += double(ap2)/double(nr);
		data[k*nfields+8] += double(ap1)*double(ap1)/double(nr);
		data[k*nfields+9] += double(ap2)*double(ap2)/double(nr);
		data[k*nfields+10] += double(sh1) / double(nr);
		data[k*nfields+11] += double(sh2)/double(nr);
		data[k*nfields+12] += double(sh1)*double(sh1)/double(nr);
		data[k*nfields+13] += double(sh2)*double(sh2)/double(nr);
		if(output || k==frames){
		    write(k+f_count);
		}
		t0=t2;
	}
	f_count+=frames;
	t=wtime()-t;
	char unit[256];
	sprintf(unit,"seconds");
	if(t>300) {
		t/= 60; sprintf(unit, "minutes");
		if(t>60) {t/=60; sprintf(unit, "hours");}
	}
	if(verbose){
		printf("# Total time %.6g %s\n"
		       "# Total bac count %d\n"
		       "# Red strain total %d\n"
		       "# Blue strain total %d\n"
               "# Calculate forces time %g\n"
               "# Integrate and divide time %g\n"
               "# Remap time %g\n"
                , t, unit, num_st1+num_st2, num_st1, num_st2,
                compt_times[0], compt_times[1], compt_times[2]);
	}
}

/** Carries out a bacterial growth simulation
 * until a target number of bacteria is reached. */
void bac_sim::solve(double dt_max, double wall_time_max, int target_tb, bool output, bool verbose) {
	double wallt=0,t0;
	// o_time is the increment in seconds to output
	double o_time=300.;
	int lchpt=50, tb=0,num_st1=0,num_st2=0,c_co=0;
	char unit[256];

	//tb = total_bacteria();
    double total_area=0;
	double ap1=0, ap2=0;
    double sh1=0, sh2=0;
    // number of cells with zero crowd factor
    int nzcf = 0;
	specific_bac_count(num_st1, num_st2, nzcf, ap1, ap2, total_area, sh1, sh2);
	tb = num_st1+num_st2;

    if(tb>=target_tb) return;

	if(verbose){
		printf("# Wall time (unit)\t Sim time\t frame\t nbact\t strain0(red)\t strain1(blue)\t active0\t active1\t sheaths1\t sheaths2\t num_threads\n");
		printf("# %8.4g seconds\t%.4g\t%d\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%d\n", wallt, time, 0, tb, num_st1, num_st2, ap1, ap2, sh1, sh2, threads);
	}

	if(output) {
		write(f_count);
		f_count++;
	}

	while(wallt<wall_time_max){
		t0=wtime();
		// Calculate threads
		threads=n_bact/spars->chunk_size+1;
		if(threads>max_threads) threads=max_threads;

		// Carry out simulation step using the regular timestep until
		// within range of the target time
		for(int i=0;i<lchpt;i++){
		    step_forward(dt_max);
		}

		specific_bac_count(num_st1, num_st2, nzcf, ap1, ap2, total_area, sh1, sh2);
		tb = num_st1+num_st2;
		t0 = wtime() - t0;
		wallt += t0;

		if(verbose && int(wallt / o_time) > c_co){
			c_co++;

			sprintf(unit,"seconds");
			double tmp_t = wallt;
			if(tmp_t>300) {
				tmp_t/= 60; sprintf(unit, "minutes");
				if(tmp_t>60) {tmp_t/=60; sprintf(unit, "hours");}
			}
            printf("# %8.4g %s\t%.4g\t%d\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%d\n", tmp_t, unit, time, f_count, tb, num_st1, num_st2, ap1, ap2, sh1, sh2, threads);
			if(output){
				write(f_count);
				f_count++;
			}
		}
	}

    sprintf(unit,"seconds");
    if(wallt>300) {
        wallt/= 60; sprintf(unit, "minutes");
        if(wallt>60) {wallt/=60; sprintf(unit, "hours");}
    }

	if(verbose){
        printf("# %8.4g %s\t%.4g\t%d\t%d\t%d\t%d\t%d\n", wallt, unit, time, f_count, tb, num_st1, num_st2, threads);
		printf("# Current sim time %.6g\n", time);
		printf("# Current time %.6g %s\n"
		       "# Total bac count %d\n"
		       "# Red strain total %d\n"
		       "# Blue strain total %d\n", wallt, unit, tb, num_st1, num_st2);
	}
	if(output){
		write(f_count);
	}
}



/** Steps the simulation forward by a given time interval, by calculating
 * forces, integrating the bacteria states, and remapping the bacteria who have
 * cross block boundaries.
 * \param[in] dt the timestep to step forward by. */
void bac_sim::step_forward(double dt, bool verbose) {
    if(dt<0.) return;

    double ttmp = wtime();
	calculate_forces(dt, verbose);
    ttmp = wtime() - ttmp;
    compt_times[0] += ttmp;

    ttmp = wtime();
	integrate_and_divide(dt, verbose);
    ttmp = wtime() - ttmp;
    compt_times[1] += ttmp;

    ttmp = wtime();
	remap(verbose);
    ttmp = wtime() - ttmp;
    compt_times[2] += ttmp;

	time+=dt;
}

/** Adds a bacterium to the simulation.
 * \param[in] (x,y) the position of the bacterium.
 * \param[in] l the length of the bacterium.
 * \param[in] theta the rotation angle of the bacterium.
 * \param[in] type the integer type of the bacterium. */
void bac_sim::put(double x,double y,double l,double theta,int type, bool t6_on, int init_sheaths) {
    /*WARNING: ONLY USE DURING INITALIZATION! */
	int nx=int((x-ax)*xsp);if(nx<0) nx=0;if(nx>=m) nx=m-1;
	int ny=int((y-ay)*ysp);if(ny<0) ny=0;if(ny>=m) ny=n-1;
	int s=nx+m*ny;
	if(jmin>ny) jmin=ny;
	if(jmax<=ny) jmax=ny+1;
	if(co[s]==mem[s]) add_region_memory(s);
	// constructor of a bacterium
	ba[s][co[s]]=bact(spars,type,n_bact++,-1, t6_on, init_sheaths, x,y,l,theta);
	ba[s][co[s]].reset_ld(randos[0]);
    co[s]++;
}


/** Calculates the total number of bacteria in the simulation.
 * \return The total bacteria. */
int bac_sim::total_bacteria() {
	int tot=*co;
	for(int i=1;i<mn;i++) tot+=co[i];
	return tot;
}

void bac_sim::specific_bac_count(int &sp1, int &sp2, int &num_zero_crdfac, double &actPerc1, double &actPerc2, double& tot_area, double &sh1, double &sh2){
    
	int tot[2]={0,0};
	int act[2]={0,0};
    int sheaths[2] = {0,0};
	tot_area = 0;
	actPerc1 = 0; actPerc2 = 0;
    sh1 = 0, sh2 = 0;
    num_zero_crdfac = 0;

	for(int i=0;i<mn;i++){
		for(int j=0;j<co[i];j++){
		    // NOTE assume the type is 0 and 1
		    tot_area += ba[i][j].area();
		    if(ba[i][j].living) {
                tot[ba[i][j].type]++;
                sheaths[ba[i][j].type] += ba[i][j].sheaths;
                if(ba[i][j].t6_on) act[ba[i][j].type]++;
                if(ba[i][j].linear_crdfac(stop_force)==0) num_zero_crdfac++;
		    }
		}
	}
	sp1 = (double) tot[0];
    if(sp1!=0.) {
        actPerc1 = (double) act[0]/sp1;
        sh1 = (double) sheaths[0]/sp1;
    }

    sp2 = (double) tot[1];
    if(sp2!=0.) {
        actPerc2 = (double) act[1]/sp2;
        sh2 = (double) sheaths[1]/sp2;
    }
}

/** Doubles the memory allocation in a simulation block.
 * \param[in] s the block to double. */
void bac_sim::add_region_memory(int s) {
	mem[s]<<=1;
	if(mem[s]>=max_region_memory) {
		fprintf(stderr,"Memory allocation exceeded in region %d\n",s);
		exit(1);
	}
	bact* nba=new bact[mem[s]];
	for(int i=0;i<co[s];i++) nba[i]=ba[s][i];
	delete [] ba[s];
	ba[s]=nba;
}

/** Time-integrates the bacteria using an Euler step, and divides those
 * that exceed a length threshold.
 * \param[in] dt the timestep to use. */
void bac_sim::integrate_and_divide(double dt, bool verbose) {
#pragma omp parallel for num_threads(threads)
	for(int s=jmin*m;s<jmax*m;s++) {
		//double wch=0.;
		for(int q=0,cco=co[s];q<cco;q++) {
			bact *b=ba[s]+q;
			// integrate all the forces and torques, dead ones can also move and spin
			b->integrate(dt);
			int BT = b->type;
			b->death_clock(dt, spars->lysis_times[BT]);
			// lengthen and divide, consume resources only if the bacteria is alive
			// and that the lag time is over
			if (b->living && time > spars->lag_times[BT]) {
                // Linear force-dependent growth regulation
                double crd_fac = b->linear_crdfac(stop_force);

                {
                    double base_rate = spars->growth_rates[BT];
                    // Growth
                    switch(spars->growth_factor) {
                        case none: break;
                        case monod: break;
                        case crowd: base_rate *= crd_fac;   break;
                        case monod_n_crd: break;
                    }
                    b->actual_grate = base_rate;

                    base_rate *= dt*b->half_body_length();
                    b->l+=base_rate;
                }

                int j=omp_get_thread_num();

                {
                // Type VI reactions
                if(b->t6_on == false){
                        double prob_on = spars->rate_on[BT] * dt;
                        b->t6_on = gsl_rand(j) < prob_on;
                } else {
                    double base_sheath_rate = spars->sheath_rates[BT];
                    switch(spars->sheath_factor) {
                        case none: break;
                        case monod: break;
                        case crowd: base_sheath_rate *= crd_fac;   break;
                        case monod_n_crd: break;
                    }
                    b->actual_srate = base_sheath_rate;
                    base_sheath_rate *= dt;
                    double frate = spars->fire_rates[BT] * dt * b->sheaths;

                    int make_one = gsl_rand(j) < base_sheath_rate;
                    b->sheaths += make_one;
                    int fire_one =  gsl_rand(j) < frate;
                    if( fire_one ) {
                        b->sheaths -= fire_one;
                        b->num_fires += fire_one;
                        // Choose a neighbor at random to fire into
                        // There also a chance to fire into the environment,
                        // without hitting any target
                        int targInd = gsl_randint(b->num_neigh + 1, j);
                        if(targInd < b->num_neigh) {
                            int tmps = b->neighbors[targInd].s;
                            int tmpq = b->neighbors[targInd].q;
                            bact &targBact = ba[tmps][tmpq];
                            if (targBact.type != b->type) {
                                // race condition may happen, but we don't care really
                                targBact.stabbed = true;
                            } else {
                                // TODO: sister cells should get same relief from recycled material
                            }
                        }
                    }
                }
                } // End sheath scope

			}// Check alive scope
		}
	}

#pragma omp parallel for num_threads(threads)
	for(int s=jmin*m;s<jmax*m;s++) {
		double local_n;
        int j=omp_get_thread_num();
		for(int q=0;q<co[s];q++) {
			bact *b=ba[s]+q;
            b->mark_stabbed_as_dead();

            if(b->divide()) {
                if(co[s]==mem[s]){
                    add_region_memory(s);
                    // Since the memory is reallocated
                    // we need to reassign the pointer
                    b = ba[s]+q;
                }
#pragma omp critical
                {
                    local_n=n_bact++;
                }
                double dtheta0=random_kick(j);
                double dtheta1=random_kick(j);

                double pl=0.5*(b->l+DEFAULT_SIM_PARS::d_b_rad);
                double cth=cos(b->theta);
                double sth=sin(b->theta);

                double randl = gsl_gaussian_rand(spars->sigma_l, j);
                double lA = b->l + randl;
                double lB = b->l - randl;

                lA = (lA - DEFAULT_SIM_PARS::d_b_rad)*0.5;
                lB = (lB - DEFAULT_SIM_PARS::d_b_rad)*0.5;
                // add a new bacterium, with (x,y) position upwind of mother
                // draw a random number out of parents total sheath
                int sh_to_daughter1 = gsl_binomial_rand(b->sheaths, j);
                int sh_to_daughter2 = b->sheaths - sh_to_daughter1;

                /**
                  This point in time in general doesn't line up with the precise time the cell need to divide.
                  But the difference should be small. To see how far we are from
                  the correct statistic on interdivision time and birth size, we
                  do a little back tracking. NOTE, this only works exactly correctly, if we have constant growth rate.
                */

                int cid = co[s];
                ba[s][cid]=bact(spars,b->type,local_n,b->id, b->t6_on, sh_to_daughter2, b->x+(pl)*cth,b->y+(pl)*sth,lA, b->theta+dtheta0);
                ba[s][cid].lb= lA;
                ba[s][cid].t_birth = 0;
                //ba[s][cid].lb= length_at_birth;
                //ba[s][cid].t_birth = current_cell_time;
                ba[s][cid].reset_ld(randos[j]);
                ba[s][cid].reset_num_fires();

                ba[s][cid].vx=b->vx;
                ba[s][cid].vy=b->vy;
                // increment the total count of bacteria in block
                ba[s][cid].omega=b->omega;
                co[s] ++;

                // change the mother (x,y) position to downwind of original
                b->x-=(pl)*cth;
                b->y-=(pl)*sth;
                b->l=lB;
                b->lb=lB;
                b->t_birth = 0;
                //b->lb=length_at_birth;
                //b->t_birth = current_cell_time;
                b->reset_ld(randos[j]);
                b->reset_num_fires();

                b->theta+=dtheta1;
                b->sheaths = sh_to_daughter1;
            }
        }
    }
}

/** Calculates and stores all of the forces on all of the bacteria. */
void bac_sim::calculate_forces(double dt, bool verbose) {

	// Calculate the damping force on each bacterium due to its velocity
#pragma omp parallel for num_threads(threads)
	for(int s=jmin*m;s<jmax*m;s++) for(int q=0;q<co[s];q++){
		ba[s][q].damping_force(spars->damp);
		// set inhibitory forces and local density to zero
		ba[s][q].reset_inhibit_forces();
        ba[s][q].reset_num_neigh();
	}

	// Add forces and torques due to contacts between bacteria.

#pragma omp parallel for num_threads(threads)
    for(int s=jmin*m;s<jmax*m;s++) {
		for(int q=0;q<co[s];q++) {
			calculate_force(s,q,dt,verbose);
		}
    }

}

/** Calculates all of the forces and torques that are on a bacterium due to
 * contacts with its neighbors.
 * \param[in] s the block that the bacterium is in.
 * \param[in] q the index within the block.
 * \param[in] dt the timestep to use
 */
void bac_sim::calculate_force(int s,int q, double dt, bool verbose) {
	int li,ui,lj,uj,ci,cj,cij,cip,cjp;
	double xx,yy,prx,pry;
	bool bcap=false, bbcap=false;
	// b is the bacterium on which we want to compute forces
	bact &b=ba[s][q];
	b.disp(xx,yy);
	const double gx=b.x-xx;
    const double gy=b.y-yy;
	const double hx=b.x+xx;
    const double hy=b.y+yy;
	int_box(b.x,b.y,fabs(xx)+2*l_thresh,fabs(yy)+2*l_thresh,li,ui,lj,uj);

	for(cj=lj;cj<=uj;cj++) for(ci=li;ci<=ui;ci++) {

		cip=ci;cjp=cj;
		prx=0;while(cip>=m) {cip-=m;prx+=bx-ax;}
		while(cip<0) {cip+=m;prx-=bx-ax;}

		pry=0;while(cjp>=n) {cjp-=n;pry+=by-ay;}
		while(cjp<0) {cjp+=n;pry-=by-ay;}

		cij=cip+m*cjp;
		double kx,ky,lx,ly,ox,oy,qx,qy,px,py,rsq,nrsq,cfx,cfy;
		for(int qq=0;qq<co[cij];qq++) {
			// bb is the bacteria around bacterium b
			bact &bb=ba[cij][qq];
			// just unique pairs of bacteria is considered
			if(bb.id<=b.id) continue;

            bool tmp_cap;
			bb.ends(kx,ky,lx,ly,prx,pry,verbose);
			// compute the closest point on bact b, with end points
            // (gx, gy) and (hx, hy), to point (kx, ky), one end of bact bb
			// The point of shortest distance to (kx, ky) on bact b, is stored in (ox, oy)
            // The counter part on this shortest line segment is stored in (px, py)
			min_distance(gx,gy,hx,hy,kx,ky,ox,oy,tmp_cap);
			rsq=dis_sq(kx,ky,ox,oy);px=kx;py=ky;
            // Update if the caps are involved in this segment
            // Since we consider (kx, ky) then bbcap is by default true
            bcap=tmp_cap; bbcap = true;

			min_distance(gx,gy,hx,hy,lx,ly,qx,qy,tmp_cap);
			nrsq=dis_sq(lx,ly,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;
                // If this distance is shorter, we keep the new points
                // (ox, oy) is the point on bact b, (px, py) that on bact bb
                ox=qx;oy=qy;px=lx;py=ly;
                // update whether the cap area is in touch
                bcap = tmp_cap; bbcap = true;
			}

			min_distance(kx,ky,lx,ly,gx,gy,qx,qy,tmp_cap);
			nrsq=dis_sq(gx,gy,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;ox=gx;oy=gy;px=qx;py=qy;
                bcap = true; bbcap = tmp_cap;
			}

			min_distance(kx,ky,lx,ly,hx,hy,qx,qy,tmp_cap);
			nrsq=dis_sq(hx,hy,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;ox=hx;oy=hy;px=qx;py=qy;
                bcap = true; bbcap = tmp_cap;
			}
			// in the end, we find the closet distance between two bacteria
			// stored in (ox, oy) (px, py), if the bacteria are not crossed,
            // and if the shortest distance is larger than (2*b_rad), then they are not touching

            // First we get the y=m1*x+b1 function for the line of bact b
            double m1 = (hy-gy)/(hx-gx), b1 = gy - m1*gx, m1inv = 1./m1;
            // Then we know the lines perpendicular to that has the slope -1/m1
            // Then we can get coefficients for lines passing through (lx, ly) and (kx, ky)
            // that are also perpendicular to bact b
            double tmp_bk = ky + m1inv * kx, tmp_bl = ly + m1inv * lx;

            // Finding the intersect between bact b line and the perpendicular line passing (kx, ky)
            double vec_kx = (tmp_bk - b1)/ (m1+m1inv), vec_ky = m1 * vec_kx + b1;
            // subtract off (kx, ky) so it makes a vector emanating from (kx, ky)
            vec_kx -= kx; vec_ky -= ky;

            // Repeat from point (lx, ly)
            double vec_lx = (tmp_bl - b1)/(m1+m1inv), vec_ly = m1 * vec_lx + b1;
            vec_lx -= lx; vec_ly -= ly;
            // If bact bb lies on one side of the line of bact b, dot 1 > 0
            double dot1 = vec_kx * vec_lx + vec_ky * vec_ly;

            // A scenario where one line segment can be consider lying one half plane divided by another
            // but the converse is not true. So here, we do the calculation yet again but switch to
            // the half plane divided by bact bb.
            double m2 = (ky-ly)/(kx-lx), b2 = ky - m2*kx, m2inv = 1./m2;
            double tmp_bg = gy + m2inv*gx, tmp_bh = hy + m2inv*hx;
            double vec_gx = (tmp_bg - b2)/(m2 + m2inv), vec_gy = m2*vec_gx + b2;
            vec_gx -= gx; vec_gy -= gy;
            double  vec_hx = (tmp_bh - b2)/(m2+m2inv), vec_hy = m2*vec_hx + b2;
            vec_hx -= hx; vec_hy -= hy;
            double dot2 = vec_gx * vec_hx + vec_gy * vec_hy;

            // If the two bacteria cross, the dot product of the distance vectors computed from above would be negative
            bool crossed = (dot1<=0 && dot2<=0);

            // If the separation is larger than 2 radii, and they are not crossing
            // Then they are definitely not in touch
            if(rsq>two_b_rad_rsq && !crossed) continue;

			// if they are touching, they have a chance of firing at each other
            b.num_neigh ++;
            b.add_neighbor(cij, qq);

            bb.num_neigh ++;
            bb.add_neighbor(s, q);

            // compute the force per unit length
            // re-use nrsq
            double contact_length = sqrt(rsq);
            // This is linear elasticity model of contact
            double local_rep_force = F_r*(1 - contact_length/two_b_rad);
            // This is Hertzian model of contact
            //double local_rep_force = F_r*pow((1 - contact_length/two_b_rad), 1.5);

            // (ox, oy) -- point on the bact b
            // (px, py) -- point on bact bb
            if(bcap) {
                // find the closest cap center to the contact point on bact b
                b.add_inhibit_forces(local_rep_force, px, py, ox, oy, contact_length);
            }
            if(bbcap){
                // find the closest cap center to the contact point on bact bb
                // since the coordinates of bb is shifted due to periodic boundaries
                // to do mechanical interactions, we have to shift them back
                double back_shifted[4] = {ox-prx, oy-pry, px-prx, py-pry};
                bb.add_inhibit_forces(local_rep_force, back_shifted[0], back_shifted[1], back_shifted[2], back_shifted[3], contact_length);
            }
            // If they are crossed, we make sure the stop them from growing
            if (crossed) {
                b.inh_force_down = 2*stop_force;
                b.inh_force_up = 2*stop_force;
                bb.inh_force_down = 2*stop_force;
                bb.inh_force_up = 2*stop_force;
            }
			// force in x direction, pointing from ox to px
			cfx=local_rep_force*(px-ox);
			// force in y direction, pointing from oy to py
			cfy=local_rep_force*(py-oy);

            if(!crossed){
                // assume that the bacteria are squeezed equally
                // so the interface between two bacteria is at the mid point
                // between (ox, oy) and (px, py)
                bb.add_force(cfx,cfy,0.5*(px+ox)-prx,0.5*(py+oy)-pry);
                // add equal, opposite force to b
                b.add_force(-cfx,-cfy,0.5*(px+ox),0.5*(py+oy));
            } else {
                // if the bacteria are crossed, we reverse the force arm to try and separate them
                double new_center_x = 2*bb.x - 0.5*(px+ox);
                double new_center_y = 2*bb.y - 0.5*(py+oy);
                bb.add_force(cfx,cfy,new_center_x-prx, new_center_y-pry);
                b.add_force(-cfx,-cfy, new_center_x, new_center_y);
            }
		}
	}
}

/** Finds a bacterium with a particular ID. This is currently very inefficient,
 * scanning over the entire data structure, but is only called infrequently.
 * \param[in] shid the ID of the particle to find.
 * \param[out] (qq,ss) the block and index of the bacterium. */
bool bac_sim::search_id(int shid,int &ss,int &qq) {
	for(ss=0;ss<mn;ss++) for(qq=0;qq<co[ss];qq++)
		if(ba[ss][qq].id==shid) return true;
	return false;
}

/** Calculates the minimum distance of a given point to a bacterium.
 * \param[in] (gx,gy) the position of the first end of the bacterium.
 * \param[in] (hx,hy) the position of the second end of the bacterium.
 * \param[in] (kx,ky) the point to consider.
 * \param[out] (ox,oy) the position on the bacterium that is closest to the
 *		       given point. */
void bac_sim::min_distance(double gx,double gy,double hx,double hy,double kx,double ky,double &ox,double &oy,bool &cap) {
    // creating two vectors, one along the body of the current bacterium (px, py)
    // the other is from the point (gx, gy) on bacterium to the point to consider (kx, ky)
	double px=hx-gx,py=hy-gy,dis=(px*(kx-gx)+py*(ky-gy))/(px*px+py*py);
	cap=false;
	if(dis<0) {dis=0; cap=true;}
	if(dis>1) {dis=1; cap=true;}
	ox=gx+dis*px;oy=gy+dis*py;
}

/** Saves the bacteria positions as a text file containing integer IDs,
 * positions, lengths, and rotations.
 * \param[in] k the integer suffix to add to the filename. */
void bac_sim::write(int k) {
	sprintf(buf,"%s/f.%05d_nr%d",filename,k,iter);
	FILE *fp=safe_fopen(buf,"w");
	fprintf(fp,"#type x y half_len theta id pid on sheaths crowd_fac num_neigh actual_growth_rate actual_sheath_rate t_birth l_birth inter_div_time num_fires\n");
	for(int s=0;s<mn;s++) {
		for(int q=0;q<co[s];q++) {
			bact &b=ba[s][q];
			int write_type = b.type + (b.sheaths>0?0:2) + (b.living?0:4);
			fprintf(fp,"%d %g %g %16.14g %g "
                        "%d %d %d %d %g "
                        "%d %g %g %g %g %g %d\n",
                        write_type, b.x, b.y, b.l, b.theta,
                        b.id, b.pid, b.t6_on, b.sheaths, b.linear_crdfac(stop_force),
                        b.num_neigh, b.actual_grate, b.actual_srate,
                        b.t_birth, b.lb, b.td, b.num_fires);
		}
	}
	fclose(fp);
}


/** Scans all of the points in each block, and remaps if they have moved from
 * one block to another. */
void bac_sim::remap(bool verbose) {

	// For each block, introduce a second counter that counts the number of
	// particles initially in this block, without including any that have
	// moved into the block
	for(int s=m*jmin;s<m*jmax;s++) gh[s]=co[s];

	// Loop over the blocks within the counter
	for(int j=jmin,ijmax=jmax;j<ijmax;j++) for(int s=j*m,i=0;i<m;i++,s++) {
		int ni,nj,l;
		double x,y;
		l=0;

		// Scan over the bacteria that started within this block
		while(l<gh[s]) {
			bact &b=ba[s][l];

			// Get the particle's initial position
			x=b.x;y=b.y;

			// Calculate which block the particle's new position is
			// within
			ni=step_int((x-ax)*xsp);
			nj=step_int((y-ay)*ysp);

			// only delete if bacterium is dead, and it's been dissolved away
			bool del = !b.living & b.dissolved;
			// See if the bacterium is within the same block,
			// if it is and it needs to be kept around (alive or dead but not yet dissolved), skip
			if(!del&&ni==i&&nj==j) l++;
			else {

				if(x_prd) {
					while(ni<0) {ni+=m;b.x+=bx-ax;}
					while(ni>=m) {ni-=m;b.x-=bx-ax;}
				}
				else{
					if(ni<0) ni=0;
					else if(ni>=m) ni=m-1;
				}
				if(y_prd) {
					while(nj<0) {nj+=n;b.y+=by-ay;}
					while(nj>=n) {nj-=n;b.y-=by-ay;}
				}
				else{
					if(nj<0) nj=0;
					else if(nj>=n) nj=n-1;
				}

				// Check whether the particle's new position is
				// within the container, and it needs to be kept around
				if(!del&&ni>=0&&ni<m&&nj>=0&&nj<n) {

					// Calculate the index the new block
					// where the particle is to be stored,
					// and allocate memory if necessary
					ni+=m*nj;
					if(co[ni]==mem[ni]) add_region_memory(ni);

					// Extend computation range if needed
					if(jmin>nj) jmin=nj;
					if(jmax<=nj) jmax=nj+1;

					// Add the particle to the new block
					ba[ni][co[ni]++]=b;
				}

				// Delete the bacterium from the original block,
				// by copying the final particle in the block's
				// memory on top of it
				ba[s][l]=ba[s][--co[s]];

				// If a bact is moved or gets deleted, then the last bact gets
				// copied into current position, then we need to check that
				// otherwise, move on
				if(co[s]+1==gh[s]) gh[s]--;
				else l++;
			}
		}
	}
}
