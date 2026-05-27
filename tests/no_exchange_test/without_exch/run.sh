#!/bin/bash

# PROVIDE THE PATH FOR THE INNER REGION EXECUTABLES
exec="/home/physastro/jdgorfinkiel/CCPForge/GitLab/UKRmol-in/build_2025_f21/bin"
execout="/home/physastro/jdgorfinkiel/CCPForge/GitLab/UKRmol-in/build_2025_f21/bin"

ulimit -s unlimited

# CREATE DIRECTORY FOR OUTPUTS

mkdir outputs

# RUN INTEGRAL CALCULATION
cp target.scatci_integrals.inp inp
$exec/scatci_integrals inp 
mv  log_file.0 outputs/target.scatci_integrals.out
ln -s moints fort.16
ln -s moints fort.17
ln -s moints fort.22

# RUN TARGET CALCULATION
$exec/congen < target.congen.singlet.A1.inp > outputs/target.congen.singlet.A1.out
$exec/scatci < target.scatci.singlet.A1.inp > outputs/target.scatci.singlet.A1.out
$exec/congen < target.congen.triplet.A1.inp > outputs/target.congen.triplet.A1.out
$exec/scatci < target.scatci.triplet.A1.inp > outputs/target.scatci.triplet.A1.out

$exec/congen < target.congen.singlet.B1.inp > outputs/target.congen.singlet.B1.out
$exec/scatci < target.scatci.singlet.B1.inp > outputs/target.scatci.singlet.B1.out
$exec/congen < target.congen.triplet.B1.inp > outputs/target.congen.triplet.B1.out
$exec/scatci < target.scatci.triplet.B1.inp > outputs/target.scatci.triplet.B1.out

$exec/congen < target.congen.singlet.B2.inp > outputs/target.congen.singlet.B2.out
$exec/scatci < target.scatci.singlet.B2.inp > outputs/target.scatci.singlet.B2.out
$exec/congen < target.congen.triplet.B2.inp > outputs/target.congen.triplet.B2.out
$exec/scatci < target.scatci.triplet.B2.inp > outputs/target.scatci.triplet.B2.out
#
$exec/denprop < target.denprop.inp > outputs/target.denprop.out
mv fort.25 fort.266

# RUN N+1 CALCULATION
$exec/congen < scattering.congen.doublet.A1.inp > outputs/scattering.congen.doublet.A1.out
$exec/scatci < scattering.scatci.doublet.A1.inp  > outputs/scattering.scatci.doublet.A1.out
$execout/outer < scattering.outer.doublet.A1.inp >  outputs/scattering.outer.doublet.A1.out
mv fort.110 eigenp.A1
mv fort.111 xsec.A1

# DELETE FORT FILE
rm fort.29
rm inp
