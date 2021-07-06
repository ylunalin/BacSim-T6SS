include ./config.mk
lflags=-L. `gsl-config --libs`
iflags=-I. `gsl-config --cflags`

# Lists of files to be built
objs=common.o sim_params.o bac_sim.o bac_sim_init.o
headers=strain.hh bact.hh common.hh sim_params.hh bac_sim.hh
src=$(patsubst %.o,%.cc,$(objs))
execs=run_sim

all:
	$(MAKE) executables

executables: $(execs)

include Makefile.dep

depend:
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(objs) $(execs) libbac.a *.dSYM

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

libbac.a: $(objs)
	rm -f libbac.a
	ar rs libbac.a $^

run_sim: run_sim.cc libbac.a $(headers)
	$(cxx) $(cflags) $(iflags) -o $@ $< -lbac $(lflags)

.PHONY: clean depend executables
