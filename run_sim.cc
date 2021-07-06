#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include "common.hh"
#include "bac_sim.hh"
#include "sim_params.hh"

int main(int argc, char** argv) {
	if(argc<2) fatal_error("Must at least provide config file.\nUsage: ./run_sim <fn.cfg>\nSee sample_sim.cfg for an example.\n", 1);

	if(strcmp(argv[1], "-h")==0){
		printf("Usage: ./run_sim file.cfg\n");
		exit(1);
	}
	sim_params spars(argv[argc-1]);

	mkdir(spars.dirname,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    int nr = spars.nr, nentries = spars.frames+1;
	int seed = spars.seed;
    const int nfields=14;
    double t=0;

    double *data = new double[nfields*nentries];
    for(int i=0;i<nfields*nentries;i++) data[i]=0.;

	bac_sim bs(&spars);
	int seed_incr = 100000*bs.max_threads;
	bool output = (nr==1);
    if(spars.output_all == 1) output = spars.output_all;
    // Not thread safe
    for(int i=0;i<nr;i++){
		bool verbose=1;
		bs.reset(i, seed+i*seed_incr);
        // Carry out the simulation
        // Initialize the simulation
		printf("# Working on iteration %d out of %d...\n", i+1, nr);
        double lt=bs.wtime();

        // Add some seed bacteria
        printf("# Initializing ...\n");
        bs.random_init();

        if(verbose) printf("# Printouts of %dth iteration.\n", i+1);
        // Evolve system to target number of bacteria
        bs.solve(output, verbose);
        // Then solve for a given time
        bs.solve(output, verbose, data);

        lt = bs.wtime()-lt;
        t += lt;
    }

	// Print out some summary
	char buf[256];
	sprintf(buf,"hours");
    t /=3600.;

    puts("#");
    puts("#");
    printf("# Average statics over %d iterations.\n", nr);
    printf("# Frame\ttime\tnbact\tstrain0(red)\tstrain1(blue)\tSTD1\tSTD2\tactive0\tactive1\tSTD1\tSTD2\tsheahts1\tsheath2\tSTD1\tSTD2\n");

    for(int i=0;i<nentries;i++){
        int ind = i*nfields;
        printf("%d\t%.4g\t%14.10g\t%14.10g\t%14.10g\t%14.10g\t%14.10g\t"
            "%14.10g\t%14.10g\t%14.10g\t%14.10g\t"
            "%14.10g\t%14.10g\t%14.10g\t%14.10g\n",
        i, data[ind], data[ind+1], data[ind+2], data[ind+3], sqrt(data[ind+4] - data[ind+2]*data[ind+2]), sqrt(data[ind+5] - data[ind+3]*data[ind+3]),\
        data[ind+6], data[ind+7], sqrt(data[ind+8] - data[ind+6]*data[ind+6]),  sqrt(data[ind+9] - data[ind+7]*data[ind+7]), \
        data[ind+10], data[ind+11], sqrt(data[ind+12] - data[ind+10]*data[ind+10]),  sqrt(data[ind+13] - data[ind+11]*data[ind+11]));
    }

    printf("# Iterations time %.6g %s\n"
            , t, buf);
	delete [] data;
    puts("# Done with everything.");
}
