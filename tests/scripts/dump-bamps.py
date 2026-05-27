#!/usr/bin/env python3
#
# Convert SWINTERF boundary amplitudes to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-bamps.py fort.21
#
# Writes all information in the binary SWINTERF output file as text to standard output.
#
# Assumes that the file uses int64 and real64 types.
#

from array  import array
from struct import unpack
from sys    import argv

# The types below need to be compatible with CONGEN INTEGERs and REALs.
isize = 4; iform = 'i'
lsize = 8; lform = 'l'
rsize = 8; rform = 'd'

# Print array of real values as table
def tabulate (frm, ncols, data):
    nrows = len(data) // ncols + 1
    for i in range(0, nrows):
        a = i * ncols
        b = min((i + 1)*ncols, len(data))
        for j in range(a, b):
            print(frm.format(data[j]), end = '')
        print()

# Read GNU or Intel Fortran unformatted record as a binary string
def read_record(f):
    begin = -1
    data = bytes()
    while begin < 0:
        begin, = unpack(iform, f.read(isize))
        data += f.read(abs(begin))
        end, = unpack(iform, f.read(isize))
        if (abs(begin) != abs(end)):
            raise Exception(f'Bookmark mismatch: {begin} vs {end}')
    return data

# Open the file for binary reading.
f = open(argv[1], 'rb')

while True:

    try:

        begin, = unpack(iform, f.read(isize))
        keyrm, nrmset, nrec, ninfo, ndata = unpack(lform * 5, f.read(lsize * 5))
        end, = unpack(iform, f.read(isize))

        begin, = unpack(iform, f.read(isize))
        header = bytes(f.read(begin)).decode('utf-8').strip()
        end, = unpack(iform, f.read(isize))

        begin, = unpack(iform, f.read(isize))
        ntarg, nvib, ndis, nchan = unpack(lform * 4, f.read(lsize * 4))
        end, = unpack(iform, f.read(isize))

        begin, = unpack(iform, f.read(isize))
        mtot, stot, gutot, ion = unpack(lform * 4, f.read(lsize * 4))
        r, rmass = unpack(rform * 2, f.read(rsize * 2))
        end, = unpack(iform, f.read(isize))

        begin, = unpack(iform, f.read(isize))
        ismax, nstat, npole, ibut = unpack(lform * 4, f.read(lsize * 4))
        rmatr, = unpack(rform, f.read(rsize))
        end, = unpack(iform, f.read(isize))

        print('"{}"'.format(header))
        print('=' * (len(header) + 2))
        print()

        print('  set index                  (nchset) : ', nrmset)
        print('  number of all records      (nrec)   : ', nrec)
        print('  number of info records     (ninfo)  : ', ninfo)
        print('  number of data records     (ndata)  : ', ndata)
        print('  number of targets          (ntarg)  : ', ntarg)
        print('  number of vib channels     (nvib)   : ', nvib)
        print('  number of dis channels     (ndis)   : ', ndis)
        print('  number of elec channels    (nchan)  : ', nchan)
        print('  irreducible representation (mtot)   : ', mtot)
        print('  total spin multiplicity    (stot)   : ', stot)
        print('  total gerade/ungerade      (gutot)  : ', gutot)
        print('  residual charge            (ion)    : ', ion)
        print('  internuclear distance      (r)      : ', r)
        print('  reduced mass               (rmass)  : ', rmass)
        print('  multipole max multiplicity (ismax)  : ', ismax)
        print('  number of inner states     (nstat)  : ', nstat)
        print('                             (npole)  : ', npole)
        print('                             (ibut)   : ', ibut)
        print('  inner region radius        (rmatr)  : ', rmatr)
        print()

        if (ismax > 0):
            a = array(rform)
            a.frombytes(read_record(f))
            print('  multipole coefficients')
            print()
            tabulate('{:25.15e}', 6, a)
            print()

        begin, = unpack(iform, f.read(isize))
        estat = array(rform); estat.fromfile(f, nstat)
        end, = unpack(iform, f.read(isize))

        print('  R-matrix poles')
        print()
        tabulate('{:25.15e}', 6, estat)
        print()

        print('  Boundary amplitudes')
        print()

        wamp = array(rform);
        wamp.frombytes(read_record(f))
        for istat in range(0, nstat):
            print('    state ', istat + 1)
            print()
            tabulate('{:25.15e}', 6, wamp[istat*nchan : (istat + 1)*nchan])
            print()

    except:

        break
