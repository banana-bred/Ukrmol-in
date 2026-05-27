This directory contains an example of how to run a calculation without exchange
between the projectile and the target electrons for H+ + H2O as a target. The
input files corresponding to exactly the same calculation, but including the 
exchange, are also provided.
 

The calculation without exchange should be run first using the run.sh script in 
the folder without_exch. Once this calculation has run, the one with exchange
can be run for comparison. The latter is run using the run.sh script in the folder
with_exch.

Note that at the moment the no exchange calculations can only be performed using SCATCI.
This option is not yet implemented in MPI-SCATCI
