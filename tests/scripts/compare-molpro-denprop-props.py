#!/usr/bin/env python3
#
# Compare properties produced by Molpro and by Denprop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 molpro-properties.py  target.molpro.out  fort.24
#
# It is assumed that the property table in Molpro output file is always introduced
# by a line containing 'State-averaged charge density matrix'.
#
# Molpro matrix elements of a given property operator are expected to be groupped together
# on consecutive lines, with [ imgvn, irelstate, jmgvn, jrelstate ] quartets of these properties
# lexicographically ordered in ascending way. If the ascending order is interrupted,
# this is interpreted as a signal that the next property is for a different spin multiplicity.
# The properties known by this script are: DMX, DMY, DMZ, QMXX, QMYY, QMZZ, QMXY, QMYZ, QMYZ.
# Only lines containing '!MCSCF trans' (or '!MCSCF expec') are parsed for the properties.
#
# The script assumes that properties for the same spin multiplicities are written by Denprop
# and by Molpro, and that Molpro writes the properties for ascending spin multiplicities
# (i.e. if Denprop contains data for singlets and triplets, then Molpro is expected to
# write properties for singlets first, then for triplets).
#
# The last column in the comparative table printed to standard output by this script contains
# difference of absolute values, rather than simple difference between the Molpro
# and Denprop properties. The reason for this is that the signs of the properties
# (except for the diagonal ones) are not physical and more or less random, depending
# on the mood of Scatci.
#
# Tested with Molpro 2010 and Denprop from UKRmol+ 2.0.
#

import math
import re
import sys

# Require two files.
if len(sys.argv) <= 2:
    print('Too few files given')
    sys.exit(1)

# Open the provided files.
molpro_file = open(sys.argv[1], 'r')
denprop_file = open(sys.argv[2], 'r')

# Check that the open succeeded.
if molpro_file == None or denprop_file == None:
    print('Invalid files given')
    sys.exit(1)

# Initialize storage lists and dicts
#  - `denprop_targets`: List of lists [ index, mgvn, spin_multiplicity ],
#                       one for each target defined in Denprop output.
#  - `denprop_props`: Dictionary of dictionaries, where the outer key is (istate,jstate)
#                     in the Denprop indexing convention, inner key is (l,m) and inner
#                     value is the property.
#  - `molpro_props`: As above, but with values read from Molpro output. Where needed,
#                    the properties are recalculated to match denprop conventions.
denprop_targets = []
denprop_props = {}
molpro_props = {}

# Read data from the Denprop property file and populate `denprop_targets`
# and `denprop_props`. Also populate a dict of lists `absolute`, where the keys
# are unique (mgvn,spin) of denprop targets and the value is a one-based list
# of ascending absolute Denprop target indices.
print('  ist   M   S   mol')
print(' ----  --  --  ----')
absolute = {}
for line in denprop_file:
    words = line.split()
    if (len(words) == 0):
        continue
    if (words[0] == '5'):
        indx = int(words[1])
        mgvn = int(words[4])
        spin = int(words[5])
        denprop_targets.append([ indx, mgvn, spin ])
        relindx = len([state[0] for state in denprop_targets if state[1] == mgvn and state[2] == spin])
        molabel = '{:4d}.{:d}'.format(relindx, mgvn + 1)
        print('{:5d}{:4d}{:4d}'.format(indx, mgvn, spin) + molabel)
        if not (mgvn,spin) in absolute:
            absolute[mgvn,spin] = []
        absolute[mgvn,spin].append(indx)
    if (line[0] == '1'):
        istate = int(words[1])
        jstate = int(words[3])
        l      = int(words[6])
        m      = int(words[7])
        prop   = float(words[8].replace('D', 'E'))
        istate, jstate = max(istate, jstate), min(istate, jstate)
        if not (istate,jstate) in denprop_props:
            denprop_props[istate,jstate] = {}
        denprop_props[istate,jstate][l,m] = prop
print()

# Get spin multiplicities present in Denprop. These will be used also in Molpro.
spins = sorted(list(set([ state[2] for state in denprop_targets ])))

# Read data from the Molpro output file.
rex = re.compile('\<([0-9]+)\.([0-9]+)\|([A-Z]+)\|([0-9]+)\.([0-9]+)\>')
for line in molpro_file:
    if 'State-averaged charge density matrix' in line:
        molpro_props = {}
        spin = -1
        last_code = ''
        last_stat = [ 0, 0, 0, 0 ]
    if '!MCSCF trans' in line or '!MCSCF expec' in line:
        words = line.split()
        prop = float(words[3])
        code = words[2]
        match = rex.match(code)
        if match:
            istate = int(match.group(1))
            imgvn = int(match.group(2)) - 1
            code = match.group(3)
            jstate = int(match.group(4))
            jmgvn = int(match.group(5)) - 1
            if code != last_code:
                spin = 0
            elif [ imgvn, istate, jmgvn, jstate ] <= last_stat:
                spin = spin + 1
            last_code = code
            last_stat = [ imgvn, istate, jmgvn, jstate ]
            istate = absolute[imgvn,spins[spin]][istate - 1]
            jstate = absolute[jmgvn,spins[spin]][jstate - 1]
            istate, jstate = max(istate, jstate), min(istate, jstate)
            if not (istate,jstate) in molpro_props:
                molpro_props[istate,jstate] = {}
            molpro_props[istate,jstate][code] = prop

# Print comparative table, converting the Molpro properties to Denprop convention
# on the fly.
print('  ist  jst    l    m    denprop     molpro    diffabs')
print(' ---- ----  ---  ---    -------     ------     ------')
disagreement = 0
for l in [ 0, 1, 2 ]:
    for m in range(-l, l+1):
        for istate, jstate in [ (istate[0], jstate[0]) for istate in denprop_targets for jstate in denprop_targets if istate[0] >= jstate[0] ]:
            uprop = 0.0
            if (istate,jstate) in denprop_props:
                if (l,m) in denprop_props[istate,jstate]:
                    uprop = denprop_props[istate,jstate][l,m]
            mprop = 0.0
            if (istate,jstate) in molpro_props:
                if l == 1 and m == -1 and 'DMY' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['DMY']
                if l == 1 and m ==  0 and 'DMZ' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['DMZ']
                if l == 1 and m == +1 and 'DMX' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['DMX']
                if l == 2 and m == -2 and 'QMXY' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['QMXY'] * 2 / math.sqrt(3)
                if l == 2 and m == -1 and 'QMYZ' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['QMYZ'] * 2 / math.sqrt(3)
                if l == 2 and m ==  0 and 'QMZZ' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['QMZZ']
                if l == 2 and m == +1 and 'QMXZ' in molpro_props[istate,jstate]:
                    mprop = molpro_props[istate,jstate]['QMXZ'] * 2 / math.sqrt(3)
                if l == 2 and m == +2 and 'QMXX' in molpro_props[istate,jstate]:
                    mprop = mprop + molpro_props[istate,jstate]['QMXX'] / math.sqrt(3)
                if l == 2 and m == +2 and 'QMYY' in molpro_props[istate,jstate]:
                    mprop = mprop - molpro_props[istate,jstate]['QMYY'] / math.sqrt(3)
            if (uprop != 0 or mprop != 0):
                diff = abs(uprop) - abs(mprop)
                disagreement = disagreement + abs(diff)
                print("{:5d}{:5d}{:5d}{:5d}{:11.6f}{:11.6f}{:11.6f}".format(istate, jstate,  l, m, uprop, mprop, diff))
print()

# Judge similarity by sum of absolute values of differences of absolute values.
print('Sum of differences: ', disagreement)
sys.exit(0 if disagreement <= 1e-5 else 1)

