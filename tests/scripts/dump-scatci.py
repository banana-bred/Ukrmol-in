#!/usr/bin/env python3
#
# Convert SCATCI output to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-scatci.py fort.25
#
# Writes all information in the binary SCATCI output file as text to standard output.
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

# Open the file for binary reading.
f = open(argv[1], 'rb')

while True:

    try:

        # Read header.
        begin, = unpack(iform, f.read(isize))
        header = bytes(f.read(32)).decode('utf-8').strip()
        end, = unpack(iform, f.read(isize))

        # Read metadata.
        begin, = unpack(iform, f.read(isize))
        nset, nrec = unpack(lform * 2, f.read(lsize * 2))
        scatci_name = bytes(f.read(120)).decode('utf-8').strip()
        nnuc, nocsf, nstat, mgvn = unpack(lform * 4, f.read(lsize * 4))
        s, sz = unpack(rform * 2, f.read(rsize * 2))
        nelt, = unpack(lform, f.read(lsize))
        e0, = unpack(rform, f.read(rsize))
        ntgsym, = unpack(lform, f.read(lsize))
        mcont_ntgs = array(lform); mcont_ntgs.fromfile(f, ntgsym*2)
        end, = unpack(iform, f.read(isize))

        print('"{}"'.format(scatci_name))
        print('=' * (len(scatci_name) + 2))
        print()

        print('  set index                  (nset)   : ', nset)
        print('  number of records          (nrec)   : ', nrec)
        print('  number of CSFs             (nocsf)  : ', nocsf)
        print('  number of states           (nstat)  : ', nstat)
        print('  irreducible representation (mgvn)   : ', mgvn)
        print('  total spin magnitude       (s)      : ', s)
        print('  total spin projection      (sz)     : ', sz)
        print('  number of electrons        (nelt)   : ', nelt)
        print('  base energy                (e0)     : ', e0)
        print('  number of target irrs.     (ntgsym) : ', ntgsym)
        print('  irrs. of continuum orbs.   (mcont)  : ', [i for i in mcont_ntgs[0::2]])
        print('  number of cont. orbs.      (notgt)  : ', [i for i in mcont_ntgs[1::2]])
        print()

        print('  Geometry:')
        for i in range(0, nnuc):
            begin, = unpack(iform, f.read(isize))
            cname = bytes(f.read(8)).decode('utf-8').strip()
            xnuc, ynuc, znuc, charge = unpack(rform * 4, f.read(rsize * 4))
            end, = unpack(iform, f.read(isize))
            print('    {:5} {:+8f} {:+8f} {:+8f} {:2f}'.format(cname, xnuc, ynuc, znuc, charge))
        print()

        begin, = unpack(iform, f.read(isize))
        iphz = array(lform); iphz.fromfile(f, nocsf)
        eig = array(rform); eig.fromfile(f, nstat)
        dg = array(rform); dg.fromfile(f, nocsf)
        end, = unpack(iform, f.read(isize))

        print('  Phases:')
        tabulate('{:5d}', 30, iphz)
        print()

        print('  Eigenenergies (+ e0):')
        tabulate('{:25.15e}', 6, [ e + e0 for e in eig ])
        print()

        print('  Diagonals:')
        tabulate('{:25.15e}', 6, dg)
        print()

        for i in range(0, nstat):
            begin, = unpack(iform, f.read(isize))
            k, = unpack(lform, f.read(lsize))
            vec = array(rform); vec.fromfile(f, nocsf)
            end, = unpack(iform, f.read(isize))

            print('  Vec #{}:'.format(i + 1))
            tabulate('{:25.15e}', 6, vec)
            print()

    except:

        break
