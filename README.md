# BacSim-T6SS
An agent-based model to simulate bacteria with type VI secretion system

Authors:
Yuexia Luna Lin
_FlexLab, Ecole Polytechnique Federale de Lausanne, Switzerland_

Chris H. Rycroft
_Department of Mathematics, University of Wisconsin-Madison_

## Software/library requirements:

- GNUMake
- GNU GCC compiler
- GSL (GNU Scientific Library)
- Perl         (For image processing only)
- Povray 3.7   (For image processing only)
- FFMPEG       (For making movies)

## Get started:
1. Download or <code> git clone </code> BacSim-T6SS repository.
2. Modify the config.mk file to configure which compiler to use and basic compiler and linker flags. The default config.mk file provide Mac Ports installed GCC 9 on Macs, another sample config for Linux systems is provided in the configs/ folder.
3. After these changes, type <code> make </code> in the command line in BacSim-T6SS/ directory, which will build an executable run_sim. This executable will use OpenMP to parallelize the computation if it is available.
4. To run the executable, a simulation config file (must have file extension .cfg) for the simulation must be provided as the only command line argument.
5. Sample simulation config files can be found in sim_configs/ directory. Examples include
 * si-movie-1.cfg : competition between wildtype ES401 and wildtype FQ-A002, where both lethal strains start from liquid culture T6SS activity level. The simulation uses periodic boundary condition to mimic the interior of a colony.
 * si-movie-2.cfg : competition between wildtype ES401 and wildtype FQ-A002, where both lethal strains start with fully activated T6SS. Boundary condition is the same as above.
 * si-movie-3.cfg : competition between a T6SS+ lethal strain and a nonlethal strain that does not fire any T6SS attacks, in a range expansion.
 
These examples produce the data files needed to render the corresponding supplementary movies in the manuscript listed below.

6. If POV-Ray is installed and the necessary output files (f.%05d_nr%d) are available, the perl script pov-movie.pl can be used to render the output files.
For example, if we want to render every frame in the directory my_sim.out/ and then link them into a movie, we can issue the command <pre> perl pov-movie.pl -m my_sim.out</pre> where -m is a flag that tells the script to call ffmpeg and link the rendered pngs into a movie.
In the pov_headers/ directory, we also provide some header files that define the colors, textures, and other rendering features for POV-Ray. For detailed usage of pov-movie.pl, see the help page by typing <pre> perl pov-movie.pl -h </pre>

More complete documentation is under development.
In the meantime, for any questions, feel free to contact the authors of this repository.

## Acknowledgement
This work has been partially supported by the Applied Mathematics Program of the U.S. DOE Office of Science Advanced Scientific Computing Research under contract number DE-AC02-05CH11231, the Department of Energy Computational Science Graduate Fellowship, and the Harvard NSF-Simons Center Quantitative Biology Initiative student fellowship, supported by the National Science Foundation grant DMS-1764269.

## References
Y.L.L, S.N.S., E.K., A.N.S, and C.H.R. "A subcellular biochemical model for T6SS dynamics reveals winning competitive strategies," in preparation.

