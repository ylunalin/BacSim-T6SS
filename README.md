# BacSim-T6SS
An agent-based model to simulate bacteria with type VI secretion system

Authors: Yuexia Luna Lin, Chris H. Rycroft

_John A. Paulson School of Engineering and Applied Sciences, Harvard University_
## Software/library requirements:

- GNUMake
- GNUGCC compiler
- GSL (GNU Scientific Library)
- Perl         (For image processing only)
- Povray 3.7   (For image processing only)

## Get started:
1. Download or <code> git clone </code> BacSim-T6SS repository.
2. Modify the config.mk file to configure which compiler to use and basic compiler and linker flags. Sample config.mk files are provided for Mac BacSim-T6SS/ folder.
3. After these changes, type <code> make </code> in the commandline in BacSim-T6SS/ directory, which will build an executable run_sim. This executable will use OpenMP to parallelize the computation if it is available.
4. To run the executable, a config file (must have file extension .cfg) for the simulation must be provided as command line argument.
5. Sample simulation config files can be found in sim_configs/ directory. Examples include
 * 
6. If POV-Ray is installed and the necessary output files (f.%05d_nr%d) are available, the perl script pov-movie.pl can be used to render the output files. For detailed usage, see the output of <pre> perl pov-movie.pl -h </pre>

More complete documentation is under development.
In the meantime, for any questions, feel free to contact the authors of this repository.

## Acknowledgement
This work has been partially supported by the Applied Mathematics Program of the U.S. DOE Office of Science Advanced Scientific Computing Research under contract number DE-AC02-05CH11231, the Department of Energy Computational Science Graduate Fellowship, and the Harvard NSF-Simons Center Quantitative Biology Initiative student fellowship, supported by the National Science Foundation grant DMS-1764269.

## References
Y.L.L, S.N.S., E.K., A.N.S, and C.H.R. "A subcellular biochemical model for T6SS dynamics reveals winning competitive strategies," in preparation.

