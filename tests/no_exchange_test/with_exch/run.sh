#!/bin/bash

# PROVIDE THE PATH FOR THE INNER REGION EXECUTABLES
exec="/home/physastro/jdgorfinkiel/CCPForge/GitLab/UKRmol-in/build_2025_f21/bin"
execout="/home/physastro/jdgorfinkiel/CCPForge/GitLab/UKRmol-in/build_2025_f21/bin"

ulimit -s unlimited

# CREATE DIRECTORY FOR OUTPUTS

mkdir outputs

# link to files created in the no exchange calculation

ln -s ../moints fort.16
ln -s ../moints fort.17
ln -s ../moints fort.22
ln -s ../fort.24 fort.24
ln -s ../fort.266 fort.266

# RUN N+1 CALCULATION
$exec/congen < scattering.congen.doublet.A1.inp > outputs/scattering.congen.doublet.A1.out
$exec/scatci < scattering.scatci.doublet.A1.inp  > outputs/scattering.scatci.doublet.A1.out
$execout/outer < scattering.outer.doublet.A1.inp >  outputs/scattering.outer.doublet.A1.out
mv fort.110 eigenp.A1
mv fort.111 xsec.A1

# DELETE FORT FILE
rm fort.29
rm inp
