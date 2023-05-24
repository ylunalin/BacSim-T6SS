#ifndef BAC_SIM_HH
#define BAC_SIM_HH

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#include "common.hh"
#include "file_output.hh"
#include "strain.hh"
#include "bact.hh"
#include "sim_params.hh"

/** The initial memory per each simulation block. */
const int init_region_memory=32;

/** The maximum memory that can be allocated to a particular simulation. */
const int max_region_memory=65536;

class bac_sim {
	public:
		/** A simulation parameter object. */
		sim_params *spars;
		/** Periodicity in x,y directions. */
		const bool x_prd, y_prd;
		/** The number of blocks in the x direction. */
		const int m;
		/** The number of blocks in the y direction. */
		const int n;
		/** The total number of blocks. */
		const int mn;
		/** The minimum coordinate in the x direction. */
		const double ax;
		/** The maximum coordinate in the x direction. */
		const double bx;
		/** The minimum coordinate in the y direction. */
		const double ay;
		/** The maximum coordinate in the y direction. */
		const double by;
		/** The length of a block in the x direction. */
		const double dx;
		/** The length of a block in the y direction. */
		const double dy;
		/** The inverse length of a block in the x direction. */
		const double xsp;
		/** The inverse length of a block in the y direction. */
		const double ysp;
		/** A global variable for crowding factor based on carrying capacity */
        double global_crd_fac;
		/** The name of the directory in which to store the output. */
		const char *filename;
		/** A counter for the total number of bacteria inserted, used
		 * to initialize the bacteria ID. (If bacteria are removed from
		 * the simulation, then n_bact may be higher than the actual
		 * number of bacteria. */
		int n_bact;
		/** A counter for the simulation frame to output. */
		int f_count;
		/** The number of iterations. */
		int nr;
		/** The simulation time. */
		double time;
		/** An array containing the number of bacteria per block. */
		int *co;
		/** An array used to keep track of the number of bacteria per
		 * block during the remapping phase. */
		int *gh;
		/** An array containing the memory allocation per block. */
		int *mem;
		/** An array container the resource field per block. */
		// double *w;
		/** An array of bacteria information per block. */
		bact **ba;

        /** Two strains */
        strain bac_types[2];
        double compt_times[3];

		/** Constructors */
		bac_sim(sim_params *sp);
		~bac_sim();

		void reset(int iter_, int seed_);
		// Time stepping function users can call
		void solve(bool output, bool verbose);
		void solve(bool output, bool verbose, double *data);
		// Under the hood solve functions
		void solve(double dt_max, double wall_time_max, int target_tb, bool output, bool verbose);
		void solve(double duration,int frames,double dt_max, bool output, bool verbose, double *data);
		// Time integration functions
		void step_forward(double dt, bool verbose=false);
		void integrate_and_divide(double dt, bool verbose=false);
		void calculate_forces(double dt,bool verbose=false);
		void remap(bool verbose=false);
		void put(double x, double y, double l, double theta, int type, bool t6_on, int init_sheaths);

		// Output and data collections
		void write(int k);

		/** Initialization functions */
		inline void random_init(int nbact, double r1, double r2, int type, double ptype0, double init_mean_sh[2]){
			if(type==0) random_circle(nbact, r1, ptype0, init_mean_sh);
			else if(type==1) random_square(nbact, r1, ptype0, init_mean_sh);
			else if(type==2) line(nbact,r1, ptype0, init_mean_sh);
			else if(type==3) half_plane(nbact, r1, ptype0, init_mean_sh);
			else if(type==4) random_square_w_circle(nbact, r1, r2, init_mean_sh);
			else {
				printf("bac_sim:: radom_init(int, double, int): Unknown init type %d\n", type);
				exit(1);
			}
		}
		inline void random_init(){
			if(internal_sp){
				printf("Warning: There is an internal sim_params object with default options."
					   " Are you sure you want to use default initial nbact, region size r, type, frames params?\n");
			}
			random_init(spars->n_bac, spars->init_region_size, spars->init_sec_lengthscale, spars->itype, spars->ptype0, spars->init_mean_sh);
		}
		void random_circle(int nbact, double r, double pty0, double init_mean_sh[2]);
		void random_square(int nbact, double r, double pty0,  double init_mean_sh[2]);
		void line(double nbact, double r, double pyt0,  double init_mean_sh[2]);
		void circle(double dh, double r, int ring_layer, double rdh, double pty0,  double init_mean_sh[2]);
		void half_plane(int nbact, double l,  double pty0, double init_mean_sh[2]);
		void random_square_w_circle(int nbact, double l, double r,  double init_mean_sh[2]);

		/** Data analysis functions */
		int total_bacteria();
        double cal_global_crd_fac();
		void specific_bac_count(int &sp1, int &sp2, int &num_zero_crd_fac, double &ap1, double &ap2, double &tot_area, double &sh1, double &sh2);

	//private:
		/** If bac_sim has an internal sim_params object */
		const bool internal_sp;
		/** The current minimum j index where a bacterium is. */
		int jmin;
		/** The current maximum j index where a bacterium is. */
		int jmax;
		/** The total number of threads to currently use. For small
		 * numbers of bacteria, maximum efficiency is achieved for a
		 * single thread. As more bacteria are introduced, this is
		 * increased. */
		int threads;
		/** A counter to differentiate between output files from
			multiple realizations. */
		int iter;
		/** A counter for diffusion solve */
		int tot_l;
		/** The maximum number of threads to use. */
		const int max_threads;
        /** Diameter of the cell, and diameter */
        const double two_b_rad;
        /** A useful constant in determining bacteria overlaps. */
        const double two_b_rad_rsq;
		/** Threshold to determine if two bac overlap.*/
		const double l_thresh;
		/** Repulsive force constant. */
        double F_r;

		/** Nondimensional growth yield */
		const double gamma;
		const double gsp;

		/** Stopping force fraction (of the repulsion force). */
		double stop_force;

		/** A temporary character buffer used to create output
		 * filenames. */
		char *buf;
		/** The thread-based random number generators. */
		gsl_rng** const randos;
		void add_region_memory(int s);
		bool search_id(int shid,int &ss,int &qq);
		void calculate_force(int s,int q, double dt,bool verbose=false);
		void min_distance(double gx,double gy,double hx,double hy,double kx,double ky,double &ox,double &oy,bool &cap);
		bool check_overlap(double xpos, double ypos, double ang, double initl);
		/** Custom int function, that gives consistent stepping for
		 * negative numbers. With normal int, we have
		 * (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1). With this routine, we
		 * have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
		inline int step_int(double a) {
			return a<0?int(a)-1:int(a);
		}
		inline void int_box(double x,double y,double r,int &li,int &ui,int &lj,int &uj) {
			int_box(x,y,r,r,li,ui,lj,uj);
		}
		inline void int_box(double x,double y,double rx,double ry,int &li,int &ui,int &lj,int &uj) {
			li=step_int((x-ax-rx)*xsp);ui=int((x-ax+rx)*xsp);
			if(!x_prd) {
				if(li<0) li=0;
				if(li>=m) li=m-1;
				if(ui<0) ui=0;
				if(ui>=m) ui=m-1;
			}
			lj=step_int((y-ay-ry)*ysp);uj=int((y-ay+ry)*ysp);
			if(!y_prd) {
				if(lj<0) lj=0;
				if(lj>=n) lj=n-1;
				if(uj<0) uj=0;
				if(uj>=n) uj=n-1;
			}
		}
		inline double random_kick(int thr_num) {
			return gsl_ran_gaussian(randos[thr_num], spars->sigma_a);
		}
		inline void set_rnd(int seed){
			for(int i=0;i<max_threads;i++) gsl_rng_set(randos[i], seed+i);
		}
		inline double dis_sq(double gx,double gy,double qx,double qy) {
			double ex=gx-qx,ey=gy-qy;
			return ex*ex+ey*ey;
		}
		inline double gsl_rand(int thr_num) const {
			return gsl_rng_uniform(randos[thr_num]);
		}
		inline double gsl_randint(int num, int thr_num) const {
			return gsl_rng_uniform_int(randos[thr_num], (unsigned long) num);
		}
		inline double gsl_r_ang(int thr_num) const{
			return gsl_rng_uniform(randos[thr_num])*2*M_PI;
		}
        inline double gsl_gaussian_rand(double sig, int thr_num) const{
            return gsl_ran_gaussian(randos[thr_num], sig);
        }
        inline double gsl_poisson_rand(double mu, int thr_num) const {
            return gsl_ran_poisson(randos[thr_num], mu);
        }
        inline double gsl_binomial_rand(int N, int thr_num) const{
            return gsl_ran_binomial(randos[thr_num], 0.5, N);
        }
        double random_length(int type, int thr_num) const{
            double rnd = gsl_rand(thr_num);
            if(spars->sigma_l ==0) rnd = 0;
            return (1+rnd)*spars->l0s[type] + rnd*r0;
        }
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return 0;}
#endif
};

#endif
