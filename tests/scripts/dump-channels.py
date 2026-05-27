#!/usr/bin/env python3
#
# Convert SWINTERF channel table to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-channels.py fort.10
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

# Open the file for binary reading.
f = open(argv[1], 'rb')

while True:

    try:

        begin, = unpack(iform, f.read(isize))
        keych, nchset, nrec, ninfo, ndata = unpack(lform * 5, f.read(lsize * 5))
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

        print('"{}"'.format(header))
        print('=' * (len(header) + 2))
        print()

        print('  set index                  (nchset) : ', nchset)
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

        print()
        print('  Target states')
        print()

        for itarg in range(0, ntarg):

            begin, = unpack(iform, f.read(isize))
            k, mtarg, starg, gutarg = unpack(lform * 4, f.read(lsize * 4))
            etarg, = unpack(rform, f.read(rsize))
            end, = unpack(iform, f.read(isize))

            print(f'    {k:4d}: {mtarg} {starg} {gutarg} {etarg:25.15e}')

        print()
        print('  Channels')
        print()

        for ichan in range(0, nchan):

            begin, = unpack(iform, f.read(isize))
            i, ichl, lchl, mchl = unpack(lform * 4, f.read(lsize * 4))
            echl, = unpack(rform, f.read(rsize))
            end, = unpack(iform, f.read(isize))

            print(f'    {i:4d}: {ichl:4d} {lchl:2d} {mchl:+3d} {echl:25.15e}')

        if (nvib > 0 or ndis > 0):

            print()
            print('  Vibrational and dissociation channels')
            print()

            for ivdch in range(0, nvib + ndis):

                begin, = unpack(iform, f.read(isize))
                i, ivtarg, iv = unpack(lform * 2, f.read(lsize * 2))
                end, = unpack(iform, f.read(isize))

                print(f'    {i:4d}: {ivtarg} {iv}')

        print()

    except:

        break
