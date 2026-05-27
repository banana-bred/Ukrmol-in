#!/usr/bin/env python3
#
# Convert CONGEN output to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-congen.py fort.70
#
# Writes all information in the binary CONGEN output file as text to standard output.
#
# Assumes that the file uses int64 and real64 types.
#

from array  import array
from struct import unpack
from sys    import argv

# Expect electron (= 0) or positron (= 1) configuration set.
iposit = 0

# The types below need to be compatible with CONGEN INTEGERs and REALs.
isize = 4; iform = 'i'
lsize = 8; lform = 'l'
rsize = 8; rform = 'd'

# Open the file for binary reading.
f = open(argv[1], 'rb')

# Read header.
begin, = unpack(iform, f.read(isize))
congen_name = bytes(f.read(120)).decode('utf-8').strip()
mgvn, = unpack(lform, f.read(lsize))
s, sz, r, pin = unpack(rform * 4, f.read(rsize * 4))
norb, nsrb, nocsf, nelt, lcdof, idiag, nsym, symtyp, lndof = unpack(lform * 9, f.read(lsize * 9))
npflag = array(lform); npflag.fromfile(f, 6)
thresh, = unpack(rform, f.read(rsize))
nctarg, ntgcon = unpack(lform * 2, f.read(lsize * 2))
end, = unpack(iform, f.read(isize))

print('"{}"'.format(congen_name))
print('=' * (len(congen_name) + 2))
print()

print('irreducible representation (mgvn)   : ', mgvn)
print('total spin magnitude       (s)      : ', s)
print('total spin projection      (sz)     : ', sz)
print('reflection symmetry        (r)      : ', r)
print('pin                        (pin)    : ', pin)
print('number of orbitals         (norb)   : ', norb)
print('number of spin-orbitals    (nsrb)   : ', nsrb)
print('number of CSFs             (nocsf)  : ', nocsf)
print('number of electrons        (nelt)   : ', nelt)
print('number of determinants     (lcdof)  : ', lcdof)
print('diagonalization flag       (idiag)  : ', idiag)
print('total number of irrs.      (nsym)   : ', nsym)
print('point group class          (symtyp) : ', symtyp)
print('length of det storage      (lndof)  : ', lcdof)
print('print flags                (npflag) : ', [flag for flag in npflag])
print('threshold                  (thresh) : ', thresh)
print('number of CI targets       (nctarg) : ', nctarg)
print('number of CI target irrs.  (ntgcon) : ', ntgcon)
print()

print('positron flag - input      (iposit) : ', iposit)
print()

if (iposit == 0):
    if (nctarg > 0):
        if (ntgcon > 0):
            begin, = unpack(iform, f.read(isize))
            iphz = array(lform); iphz.fromfile(f, nctarg)
            nctgt = array(lform); nctgt.fromfile(f, ntgcon)
            notgt = array(lform); notgt.fromfile(f, ntgcon)
            mcont = array(lform); mcont.fromfile(f, ntgcon)
            gucont = array(lform); gucont.fromfile(f, ntgcon)
            end, = unpack(iform, f.read(isize))
            print('taget phase correction     (iphz)   : ', [phz for phz in iphz])
            print('number of targets per irr  (nctgt)  : ', [n for n in nctgt])
            print('number of cont. orbs.      (notgt)  : ', [n for n in notgt])
            print('symmetry of cont. orbs.    (mcont)  : ', [m for m in mcont])
            print('reflection sym. of orbs.   (gucont) : ', [gu for gu in gucont])
            print()
        else:
            begin, = unpack(iform, f.read(isize))
            iphz = array(lform); iphz.fromfile(f, nctarg)
            print('taget phase correction     (iphz)   : ', [phz for phz in iphz])
            print()
            end, = unpack(iform, f.read(isize))

        begin, = unpack(iform, f.read(isize))
        nob = array(lform); nob.fromfile(f, nsym)
        ndtrf = array(lform); ndtrf.fromfile(f, nelt)
        nodo = array(lform); nodo.fromfile(f, nocsf)
        iposit, = unpack(lform, f.read(lsize))
        nob0 = array(lform); nob0.fromfile(f, nsym)
        nobl = array(lform); nobl.fromfile(f, 2*nsym)
        nob0l = array(lform); nob0l.fromfile(f, 2*nsym)
        end, = unpack(iform, f.read(isize))
        print('number of orbs. per irr.   (nob)    : ', [n for n in nob])
        print('reference determinant      (ndtrf)  : ', [n for n in ndtrf])
        print('number of dets. per CSF    (nodo)   : ', [n for n in nodo])
        print('positron flag              (iposit) : ', iposit)
        print('number of target orbitals  (nob0)   : ', [n for n in nob0])
        print('Dinfh version of nob       (nobl)   : ', [n for n in nobl])
        print('Dinfh version of nob0      (nob0l)  : ', [n for n in nob0l])
        print()
else:
    pass # TBD, see scatci_routines::rdnfto

begin, = unpack(iform, f.read(isize))
icdo = array(lform); icdo.fromfile(f, nocsf + 1)
indo = array(lform); indo.fromfile(f, nocsf + 1)
end, = unpack(iform, f.read(isize))
begin, = unpack(iform, f.read(isize))
ndo = array(lform); ndo.fromfile(f, lndof)
end, = unpack(iform, f.read(isize))
begin, = unpack(iform, f.read(isize))
cdo = array(rform); cdo.fromfile(f, lcdof)
end, = unpack(iform, f.read(isize))

for icsf in range(0, nocsf):
    print('CSF #{}'.format(icsf + 1))
    i = indo[icsf]-1
    for idet in range(icdo[icsf], icdo[icsf+1]):
        nr = ndo[i]
        src = [j for j in ndo[i+1:i+nr+1]]
        dst = [j for j in ndo[i+nr+1:i+2*nr+1]]
        print('    {:+8f} det #{:03d}: {} -> {}'.format(cdo[idet-1], idet, src, dst))
        i = i + 2*nr + 1
    print()
