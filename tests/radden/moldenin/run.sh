#!/bin/bash

# PROVIDE THE PATH FOR THE INNER REGION EXECUTABLES
exec="../../UKRmol-in/bin/"

cat > input << EOF
&INPUT
rstart=0d0        !radial interval in a.u. for which the integration will be carried out
rfinish=20d0
whichdm=1          !1=Molden file will be processed by default
iprint=.true.
&END
&MOLDENIN
DMUNIT=451                                !default input/output unit for the density matrices
molden_file='./h2o.smallest.molden'                    !path to the Molden file
header_base='h2o test'                    !default header of the density matrices if they are written to the DMUNIT, otherwise not used
&END
EOF

$exec/radden < input > output
rm ./input

