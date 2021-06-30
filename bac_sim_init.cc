#include "bac_sim.hh"


/** Check if a new bacterium overlaps with an existing one.
 * Taking into account periodic boundary conditions.
 */
bool bac_sim::check_overlap(double xpos, double ypos, double ang, double initl){
	bact stick (spars, 0, 0, 0, true, 0, xpos,ypos, initl, ang);
	int li, ui, lj, uj;
	int ci, cj, cij, cip, cjp;
	double gx,gy,hx,hy,xx,yy;
	stick.disp(xx,yy);

	gx = stick.x-xx; gy = stick.y-yy;
	hx = stick.x+xx; hy = stick.y+yy;

    double thresh_sq = two_b_rad_rsq;
	int_box(xpos, ypos, fabs(xx)+2*l_thresh, fabs(yy)+2*l_thresh,li,ui,lj,uj);

	double kx, ky, lx,ly,prx,pry;
	double ox, oy, qx, qy;
	double rsq, nrsq;

	bool loverlap=false;
    bool cap;
	for(ci=li;ci<=ui;ci++) for(cj=lj;cj<=uj;cj++){
		cip =ci; cjp = cj;
		prx=0; pry=0;
		while(cip>=m){cip-=m; prx+=bx-ax;}
		while(cip<0){cip+=m; prx-=bx-ax;}

		while(cjp>=n) {cjp-=n; pry+=by-ay;}
		while(cjp<0){cjp+=n; pry-=by-ay;}

		cij = cip + m*cjp;

		for(int i=0;i<co[cij];i++){
			bact &b=ba[cij][i];
			b.ends(kx,ky,lx,ly,prx,pry);

			min_distance(gx,gy,hx,hy,kx,ky,ox,oy,cap);
			rsq=dis_sq(kx,ky,ox,oy);

			min_distance(gx,gy,hx,hy,lx,ly,qx,qy,cap);
			nrsq=dis_sq(lx,ly,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;
			}
	
			min_distance(kx,ky,lx,ly,gx,gy,qx,qy,cap);
			nrsq=dis_sq(gx,gy,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;
			}

			min_distance(kx,ky,lx,ly,hx,hy,qx,qy,cap);
			nrsq=dis_sq(hx,hy,qx,qy);
			if(nrsq<rsq) {
				rsq=nrsq;
			}

			if(rsq<thresh_sq) {
				loverlap = true; break;
			}


            // check crossing of two bacteria
			b.ends(kx,ky,lx,ly,prx,pry);
            stick.ends(gx, gy, hx, hy, 0, 0);

            double m1 = (hy-gy)/(hx-gx), b1 = gy - m1*gx, m1inv = 1./m1;
            double tmp_bk = ky + m1inv * kx, tmp_bl = ly + m1inv * lx;
            double vec_kx = (tmp_bk - b1)/ (m1+m1inv), vec_ky = m1 * vec_kx + b1;
            vec_kx -= kx; vec_ky -= ky;

            double vec_lx = (tmp_bl - b1)/(m1+m1inv), vec_ly = m1 * vec_lx + b1;
            vec_lx -= lx; vec_ly -= ly;
            double dot1 = vec_kx * vec_lx + vec_ky * vec_ly;

            double m2 = (ky-ly)/(kx-lx), b2 = ky - m2*kx, m2inv = 1./m2;
            double tmp_bg = gy + m2inv*gx, tmp_bh = hy + m2inv*hx;
            double vec_gx = (tmp_bg - b2)/(m2 + m2inv), vec_gy = m2*vec_gx + b2;
            vec_gx -= gx; vec_gy -= gy;
            double  vec_hx = (tmp_bh - b2)/(m2+m2inv), vec_hy = m2*vec_hx + b2;
            vec_hx -= hx; vec_hy -= hy;
            double dot2 = vec_gx * vec_hx + vec_gy * vec_hy;

            if(dot1 <0 && dot2 <0) {
                loverlap = true;
                //printf("found a cross!\n%g %g  %g %g\n%g %g  %g %g\n", gx, gy, kx, ky, hx, hy, lx, ly);
                //printf("%g %g %g %g\n", vec_kx+kx, vec_ky+ky, vec_gx+gx, vec_gy+gy);
                //printf("%g %g %g %g\n", vec_lx+lx, vec_ly+ly, vec_hx+hx, vec_hy+hy);
                break;
            }
		}
	}
	return loverlap;
}

/** Intialize a random array of bacteria */
void bac_sim::random_circle(int nbact, double r, double pty0, double init_mean_sh[2]){
	int c = 0;
	double diam = r*2;
	double rsq = r*r;
    double t0 = wtime();
	while(c<nbact){
		double xpos = r;
		double ypos = r;
		double ang = 0;
        int type = c >= int (nbact*pty0);
        double initl = random_length(type, 0);
        int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
        bool on = (gsl_rand(0) < spars->init_prob_on[type]);
        if (!on) { init_sheaths = 0;}
		bool overlap = true;
		while(overlap || xpos*xpos + ypos*ypos > rsq){
			xpos = gsl_rand(0)*diam - r;
			ypos = gsl_rand(0)*diam - r;
			ang = gsl_r_ang(0);
			overlap = check_overlap(xpos, ypos, ang, initl);
            if(wtime() - t0 > 600.) {
                printf("Initialization takes too long, you sure the initialization domain can accomodate the number of bacteria you asked for?\n");
                exit(1);
            }
		}
		put(xpos,ypos,initl,ang,type,on,init_sheaths);
		c++;
	}

}

void bac_sim::half_plane(int nbact, double l, double pty0, double init_mean_sh[2]){
	int c=0;
    double t0 = wtime();
	while(c<nbact){
		double xpos=0, ypos=0, ang=0;
        int type = c >= int (nbact*pty0);
        double initl = random_length(type, 0);
        // HACK HACK HAC
        initl = spars->l0s[0];
        int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
        bool on = (gsl_rand(0) < spars->init_prob_on[type]);
        if (!on) { init_sheaths = 0;}
		bool overlap =true;
		while(true){
			xpos = gsl_rand(0)*2*l-l;
            if(type) {
                ypos = gsl_rand(0)*(-l);
            } else {
                ypos = gsl_rand(0)*l;
            }
			ang = gsl_r_ang(0);
			overlap = check_overlap(xpos, ypos, ang, initl);
			if(!overlap) break;
            if(wtime() - t0 > 600.) {
                printf("Initialization takes too long, you sure the initialization domain can accomodate the number of bacteria you asked for?\n");
                exit(1);
            }
		}
        put(xpos,ypos,initl,ang,type,on,init_sheaths);
		c++;
	}

}

void bac_sim::random_square_w_circle(int nbact, double l, double r, double init_mean_sh[2]){
	int c=0;
    double t0 = wtime();
	while(c<nbact){
		double xpos=0, ypos=0, ang=0;
        int type = sqrt(xpos*xpos + ypos*ypos)<r?0:1;
        double initl = random_length(type, 0);
        int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
        bool on = (gsl_rand(0) < spars->init_prob_on[type]);
        if (!on) { init_sheaths = 0;}
		bool overlap =true;
		while(true){
			xpos = gsl_rand(0)*2*l-l;
			ypos = gsl_rand(0)*2*l-l;
			ang = gsl_r_ang(0);
			overlap = check_overlap(xpos, ypos, ang, initl);
			if(!overlap) break;
            if(wtime() - t0 > 600.) {
                printf("Initialization takes too long, you sure the initialization domain can accomodate the number of bacteria you asked for?\n");
                exit(1);
            }
		}
        put(xpos,ypos,initl,ang,type,on,init_sheaths);
		c++;
	}
}

/** Initialize a random array of bacteria in square region */
void bac_sim::random_square(int nbact, double l, double pty0, double init_mean_sh[2]){
	int c=0;
    double t0 = wtime();
	while(c<nbact){
		double xpos=0, ypos=0, ang=0;
        int type = c >= int (nbact*pty0);
        double initl = random_length(type, 0);
        int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
        bool on = (gsl_rand(0) < spars->init_prob_on[type]);
        if (!on) { init_sheaths = 0;}
		bool overlap =true;
		while(true){
			xpos = gsl_rand(0)*2*l-l;
			ypos = gsl_rand(0)*2*l-l;
			ang = gsl_r_ang(0);
            if(nbact==1) {
                xpos = 0; ypos = 0;
            }
            else if(nbact==2) {
                initl = spars->l0s[0];
                if(c==0) {
                xpos = -l/2; ypos = -l/2;
                } else {
                xpos = l/2; ypos = l/2;
                }
            }
			overlap = check_overlap(xpos, ypos, ang, initl);
			if(!overlap) break;
            if(wtime() - t0 > 600.) {
                printf("Initialization takes too long, you sure the initialization domain can accomodate the number of bacteria you asked for?\n");
                exit(1);
            }
		}
        put(xpos,ypos,initl,ang,type,on,init_sheaths);
		c++;
	}
}

/** Initialize the bacteria colony as a circle on cartesian grid
 * \param[in] (dx,dy) horizontal and vertical separation
 * \param[in] r the radius (in number of bacteria)*/
void bac_sim::circle(double dh, double r, int ring_layer, double rdh, double pty0,double init_mean_sh[2]){
	int grid_ratio = (int) (dh/rdh);
	int rmax = (int) r+1;
	double din = r*dh;
	for (int i =-rmax; i<= rmax; i++){
		for(int j=-rmax;j<=rmax;j++){
			double x = (double) i;
			double y = (double) j;
            int type = int (gsl_rand(0)>pty0);
			double myr = x*x+y*y;
			if (myr<=r*r) {
				double xdist = x*dh;
				double ydist = y*dh;
                int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
                double initl = random_length(type, 0);
                bool on = (gsl_rand(0) < spars->init_prob_on[type]);
		if (!on) { init_sheaths = 0;}
                while(initl<0) {
                    initl = random_length(type, 0);
                }
				put(xdist, ydist, initl, gsl_r_ang(0),type, on, init_sheaths);
			}
		}
	}

	// outer layer, as if we have a finner grid
	rmax = (int) (grid_ratio*r+ring_layer+1);
	double dout = din + ring_layer*rdh;

	for (int i =-rmax; i<= rmax; i++){
		for(int j=-rmax;j<= rmax;j++){
			double x = (double) i*rdh;
			double y = (double) j*rdh;
			double myr = x*x + y*y;
			if(myr > din*din){
                int type = int (gsl_rand(0)>pty0);
				if (myr<=dout*dout) {
                    int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
                    double initl = random_length(type, 0);
                    bool on = (gsl_rand(0) < spars->init_prob_on[type]);
		if (!on) { init_sheaths = 0;}
                    while(initl<0) {
                        initl = random_length(type, 0);
                    }
					put(x, y, initl, gsl_r_ang(0),type,on,init_sheaths);
				}
			}
		}
	}
}


/** Initialize a line of bacteria.
 */
void bac_sim::line(double nbact, double r, double pty0, double init_mean_sh[2]){
	double dh = 2*r / nbact;
	if(dh<1.1*l_thresh) dh = 1.1*l_thresh;
	double x=-r;
	while(x<=r){
        int type = int (gsl_rand(0)>pty0);
        int init_sheaths = gsl_poisson_rand(init_mean_sh[type], 0);
        double initl = random_length(type, 0);
        bool on = (gsl_rand(0) < spars->init_prob_on[type]);
	if (!on) { init_sheaths = 0;}
        while(initl<0) {
            initl = random_length(type, 0);
        }
		put(x,0,initl,gsl_r_ang(0),type, on,init_sheaths);
		x+=dh;
	}
}
