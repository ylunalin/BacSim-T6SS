#ifndef BACT_HH
#define BACT_HH

#include "common.hh"
#include "strain.hh"
#include "sim_params.hh"
#include <vector>

struct bactCoord {
    int s; // block
    int q; // index within block
    bactCoord(int s_, int q_): s(s_), q(q_) {}
};


struct bact : public strain {
	/** The bacterium's ID. */
	int id;
	/** The bacterium's parent ID. */
	int pid;
    /** Number of sheaths. */
    int sheaths;
    /** Number of neighbors. */
    int num_neigh;
    /** Number of firing events. */
    int num_fires;

    bool t6_on;
	/** Status of the bacterium. */
	bool living;
	/** Whehter it has been stabbed this time step */
	bool stabbed;
	/** Status of the dead carcass of bacterium. */
	bool dissolved;

	/** The x position of the bacterium. */
	double x;
	/** The y position of the bacterium. */
	double y;
	/** The x velocity of the bacterium. */
	double vx;
	/** The y velocity of the bacterium. */
	double vy;
	/** The x force on the bacterium. */
	double fx;
	/** The y force on the bacterium. */
	double fy;
    /** The half length at birth. */
    double lb;
	/** The half length of the bacterium. */
	double l;
    /** The full length at division. */
    double ld;
    /** The interdivision time */
    double td;
	/** The rotation of the bacterium. */
	double theta;
	/** The angular velocity of the bacterium. */
	double omega;
	/** The torque on the bacterium. */
	double torque;

    /** The time since the bacterium is born */
    double t_birth;
    /** The time since the bacterium is killed */
    double t_death;
    /** Growth inhibiting force. */
    double inh_force_down;
    double inh_force_up;

    /** Actual growth rate */
    double actual_grate;
    /** Actual sheath rate */
    double actual_srate;
    std::vector<bactCoord> neighbors;

	bact():
        strain(), id(-1), pid(-1), sheaths(0), num_neigh(0), num_fires(0),
        living(0), stabbed(0), dissolved(0),
        x(0), y(0), vx(0), vy(0), fx(0), fy(0), lb(0), l(0), ld(0), td(0),
        theta(0), omega(0), torque(0),
        t_birth(0.), t_death(0),
        inh_force_down(0), inh_force_up(0),
        actual_grate(0), actual_srate(0)
        {}

	bact(sim_params * spars, int type_,int id_,int pid_, bool on, int sheaths_, double x_,double y_,double l_,double theta_) :
		strain(type_, spars), id(id_), pid(pid_), sheaths(sheaths_), num_neigh(0), num_fires(0),
        t6_on(on), living(true), stabbed(false), dissolved(false),
        x(x_), y(y_), vx(0), vy(0), fx(0), fy(0), lb(l_), l(l_), ld(0), td(0),
        theta(theta_), omega(0), torque(0),
        t_birth(0.0), t_death(0.0),
        inh_force_down(0), inh_force_up(0),
        actual_grate(0), actual_srate(0)
        {}

	inline void disp(double &xx,double &yy) const {
		xx=l*cos(theta),yy=l*sin(theta);
	}
	inline void ends(double &gx,double &gy,double &hx,double &hy,double prx,double pry, bool verbose=false) const {
		double xx=l*cos(theta),yy=l*sin(theta);
		gx=x+prx-xx;gy=y+pry-yy;
		hx=x+prx+xx;hy=y+pry+yy;
	}
	void body_vector (double &xx,double &yy) const {
		xx=(l)*cos(theta);
        yy=(l)*sin(theta);
	}
	void integrate(double dt) {
        t_birth+=dt;
		x+=dt*vx;
		y+=dt*vy;
		double dtt=dt/mass();
		vx+=dtt*fx;
		vy+=dtt*fy;
		theta+=dt*omega;
		omega+=dt*torque/m_inertia();
	}
    // Even though bacteria are 2D, we assume the mass is computed as a 3D spherocylinder
	double mass() const{
		return cap_area*(2*l + 1.3333333*r0);
	}
    // Even though bacteria are 2D, we assume the moment of inertia
    // is computed as a 3D spherocylinder
	double m_inertia() const{
		// This formula isn't quite right, but it probably doesn't
		// matter since there are many things that could influence the
		// rotational behavior.
        double h_rod_l = (l - r0);
        if(h_rod_l<0) h_rod_l=0;
        double rod_part = body_area_multiplier*h_rod_l*h_rod_l*h_rod_l*0.6666667;
        double cap_part = 2* cap_area * (0.4 + (h_rod_l+r0)*(h_rod_l+r0));
		return rod_part + cap_part;
	}
	void damping_force(double damp, bool verbose=false) {
		fx=-damp*vx;
		fy=-damp*vy;
		torque=-damp*omega*area();
	}
	void add_force(double cfx,double cfy,double cx,double cy,bool verbose=false) {
		double tque=(cx-x)*cfy-(cy-y)*cfx;
#pragma omp atomic
		fx+=cfx;
#pragma omp atomic
		fy+=cfy;
#pragma omp atomic
		torque+=tque;
	}
	void death_clock(double dt, double t_die){
		// if the bacterium is dead, start the count of death clock
		if (living==false){
			t_death += dt;
		}
		// if the bacterium has been dead for long enough, mark it dissolved
		if (t_death>t_die) dissolved = true;
	}
	void print_state() const{
		printf("Bact:: print_state() : ID %d pos(%g %g) forces (%g %g) torque %g velocities (%g %g) omega %g theta %g\n", id, x, y, fx, fy, torque, vx, vy, omega, theta);
	}
    inline void reset_inhibit_forces(){
        inh_force_down=0;
        inh_force_up=0;
     }
    void add_inhibit_forces(double f, double ox, double oy, double cx, double cy, double length){
        // (ox,oy) the coordinate of the contact point, on the other bacterium
        // (cx,cy) the coordinate of center of the cap that's under contact
        // (cx-x, cy-y) is the vector pointing from center of cell to cell of cap, l is its magnitdue
        // (ox-cx, oy-cy) is the vector between two contact points
        // length is the magnitude of (ox-cx,oy-cy)
        // f is the magnitude of contact force

        // Getting the projection of the force vector along the longitudinal direction of the bact
        double impact_x = cx-x, impact_y = cy-y;
        double body_x, body_y;
        body_vector(body_x, body_y);

        // Get dot product to see which way the contact force is pointing
        double dot_prod = body_x*impact_x + body_y*impact_y;
        double wt = f*fabs((impact_x)*(ox-cx) + (impact_y)*(oy-cy))/length/l;
        if(dot_prod<0) {
#pragma omp atomic
            inh_force_down += wt;
        } else {
#pragma omp atomic
            inh_force_up += wt;
        }
     }
    double linear_crdfac (double fmax) const {
        double denom = inh_force_up + inh_force_down;
        if(denom <1e-12) return 1.;
        double num = inh_force_up*inh_force_down;
        double tmp = 1-(num/denom/fmax);
        return tmp<0?0:tmp;
    }
    double old_lin_crdfac(double fmax) const {
        double num = inh_force_up+inh_force_down;
        double tmp = 1-(num/fmax);
        return tmp<0?0:tmp;
    }

    double monod_fac() const {
        return 1;
    }
    double area() const{
        return body_area_multiplier*l + cap_area;
    }
    double half_body_length() const {
        return r0 + l;
    }
    double body_length() const {
        return 2*(r0 + l);
    }
    bool divide() const {
        return (living && body_length()>=ld);
    }
    void reset_ld(const gsl_rng * rand) {
        ld = targ_length_w_noise(lb, td, rand);
    }
    void add_neighbor(int s, int q) {
        bactCoord tmp(s,q);
#pragma omp critical
{
        neighbors.push_back(tmp);
}
    }
    void reset_num_neigh() {
        num_neigh = 0;
        neighbors.clear();
    }
    void reset_num_fires() {
        num_fires = 0;
    }
    void mark_stabbed_as_dead (){
        if(stabbed) {
            living = false;
        }
    }
};

#endif
