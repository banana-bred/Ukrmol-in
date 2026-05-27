#!/usr/bin/env python3
#
# Convert `molecular_data` to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-molecular-data.py [N]
#
# Writes complete content of the `molecular_data` file (input for RMT) to standard output
# for debugging. No arguments are needed. The file `molecular_data` is expected to be
# in the working directory. If an integer number is provided on the command line, it will
# be used to limit number of rows printed in tables.
#
# Assumes that the file uses int32 and real64 types.
#

from array  import array
from struct import unpack
from sys    import argv

# The types below need to be compatible with RMT INTEGERs and REALs.
isize = 4; iform = 'i'
rsize = 8; rform = 'd'

# Limit on number of rows in tables
maxrows = int(argv[1]) if len(argv) > 1 else -1

# Print array of real values as table
def tabulate (frm, ncols, data):
    nrows = len(data) // ncols + 1
    if maxrows > 0:
        nrows = min(maxrows, nrows)
    for i in range(0, nrows):
        a = i*ncols
        b = min((i + 1)*ncols, len(data))
        for j in range(a, b):
            print(frm.format(data[j]), end = '')
        print()

# Open the file for binary reading.
f = open('molecular_data', 'rb')

# Read inner dipole moments.
for m in [ -1, 0, 1 ]:

    s, s1, s2 = unpack(iform * 3, f.read(3 * isize))
    print('[m = {}] s = {}, s1 = {}, s2 = {}'.format(m, s, s1, s2))

    iidip = array(iform); iidip.fromfile(f, s)
    print('iidip = ', [ i for i in iidip ])

    ifdip = array(iform); ifdip.fromfile(f, s)
    print('ifdip = ', [ i for i in ifdip ])

    # Cowardly skip dipstos larger than billion elements
    if (s * s1 * s2 <= 1000000000):
        dipsto = array(rform); dipsto.fromfile(f, s * s1 * s2)
        print('dipsto:')
        tabulate('{:25.15e}', 8, dipsto)
    else:
        print('dipsto skipped as it is HUGE ({} elements)\n'.format(s * s1 * s2))
        f.seek(s * s1 * s2 * rsize, 1)

ntarg, = unpack(iform, f.read(isize))
print('ntarg = {}'.format(ntarg))

# Read target coupling coefficients
for m in [ -1, 0, 1 ]:

    crlv = array(rform); crlv.fromfile(f, ntarg * ntarg)
    print('crlv:')
    tabulate('{:25.15e}', 8, crlv)

n_rg, = unpack(iform, f.read(isize))
print('n_rg = {}'.format(n_rg))

rg = array(rform); rg.fromfile(f, n_rg)
print('rg:')
tabulate('{:25.15e}', 8, rg)

lm_rg = array(iform); lm_rg.fromfile(f, 6 * n_rg)
print('lm_rg:')
tabulate('{:5d}', 40, lm_rg)

nelc, nz, lrang2, lamax, ntarg, inast, nchmx, nstmx, lmaxp1 = unpack(iform * 9, f.read(9 * isize))
print('nelc = {}, nz = {}, lrang2 = {}, lamax = {}, ntarg = {}, inast = {}, nchmx = {}, nstmx = {}, lmaxp1 = {}'
      .format(nelc, nz, lrang2, lamax, ntarg, inast, nchmx, nstmx, lmaxp1))

rmatr, bbloch = unpack(rform * 2, f.read(2 * rsize))
print('rmatr = {}, bbloch = {}'.format(rmatr, bbloch))

etarg = array(rform); etarg.fromfile(f, ntarg)
print('etarg = ', [ e for e in etarg ])
ltarg = array(iform); ltarg.fromfile(f, ntarg)
print('ltarg = ', [ l for l in ltarg ])
starg = array(iform); starg.fromfile(f, ntarg)
print('starg = ', [ s for s in starg ])

nfdm, = unpack(iform, f.read(isize))
delta_r, = unpack(rform, f.read(rsize))
print('nfdm = {}. delta_r = {}'.format(nfdm, delta_r))

r_points = array(rform); r_points.fromfile(f, nfdm + 1)
print('r_points = ', [ r for r in r_points ])

# Read per-symmetry data
for i in range(0, inast):

    lrgl, nspn, npty, nchan, mnp1 = unpack(iform * 5, f.read(5 * isize))
    print('[sym = {}] lrgl = {}, nspn = {}, npty = {}, nchan = {}, mnp1 = {}'.format(i, lrgl, nspn, npty, nchan, mnp1))

    nconat = array(iform); nconat.fromfile(f, ntarg)
    print('nconat = ', [ n for n in nconat ])

    l2p = array(iform); l2p.fromfile(f, nchmx)
    print('l2p = ', [ l for l in l2p[0:nchan] ])

    m2p = array(iform); m2p.fromfile(f, nchmx)
    print('m2p = ', [ m for m in m2p[0:nchan] ])

    eig = array(rform); eig.fromfile(f, nstmx)
    print('eig:')
    tabulate('{:25.15e}', 8, [ e for e in eig[0:mnp1] ])

    wmat = array(rform); wmat.fromfile(f, nchmx * nstmx)
    print('wmat:')
    tabulate('{:25.15e}', 8, wmat)

    cf = array(rform); cf.fromfile(f, nchmx * nchmx * lamax)
    print('cf:')
    tabulate('{:25.15e}', 8, cf)

    ichl = array(iform); ichl.fromfile(f, nchmx)
    print('ichl:')
    tabulate('{:5d}', 40, [ i for i in ichl[0:nchan] ])

    s1, s2 = unpack(iform * 2, f.read(2 * isize))
    print('s1 = {}, s2 = {}'.format(s1, s2))

    for j in range (0, nfdm):
        wmat2 = array(rform); wmat2.fromfile(f, s1 * s2)
        print('wmat2 ir = {}'.format(j))
        tabulate('{:25.15e}', 8, wmat2)
