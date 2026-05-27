#!/usr/bin/env python3
#
# Convert SCATCI hamiltonian to text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-scatci-hamiltonian.py fort.80
#
# Writes all information in the binary SCATCI hamiltonian file as text to
# standard output.
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

def print_matrix(indices, mat_elems, no_elems, no_col=5):
    for row in range(no_elems//no_col + 1):
        string = ""
        for j in range(5):
            idx1 = 2*(row*no_col+j)
            idx2 = row*no_col+j
            if idx2 == no_elems:
                break
            string += f'{indices[idx1]:5d}{indices[idx1+1]:5d}{mat_elems[idx2]:16.8e}'
        if string != "":
            print(string)
    return

elem_count = 0
end_of_file = False

f = open(argv[1], 'rb')

while not end_of_file:

    try:

        # Read header.
        begin, = unpack(iform, f.read(isize))
        nocsf, lembf, = unpack(lform*2, f.read(lsize*2))
        _, nocsf2 = unpack(lform*2, f.read(lsize*2))
        _, nosym = unpack(lform*2, f.read(lsize*2))
        _ = unpack(lform*4, f.read(lsize*4))
        nnuc, _ = unpack(lform*2, f.read(lsize*2))
        header = bytes(f.read(120)).decode('utf-8').strip()
        nhe = unpack(lform*20, f.read(lsize*20))
        dtnuc = unpack(rform*41, f.read(rsize*41))
        end, = unpack(iform, f.read(isize))

        print('"{}"'.format(header))
        print('=' * (len(header) + 2))
        print()

        print('nocsf  : ', nocsf)
        print('lembf  : ', lembf)
        print('nocsf2 : ', nocsf2)
        print('nosym  : ', nosym)
        print('nnuc   : ', nnuc)
        print('nhe    : ', nhe)
        print('dtnuc  : ', dtnuc)
        print()

        print('Hamiltonian matrix elements:')
        print()

        # Read the hamiltonian matrix elements in chunks.
        while not end_of_file:

            try:
                begin, = unpack(iform, f.read(isize))
                no_elem, = unpack(lform, f.read(lsize))
                elem_count += no_elem
                indices = unpack(lform*2*lembf, f.read(lsize*2*lembf))
                elements = unpack(rform*lembf, f.read(rsize*lembf))
                end, = unpack(iform, f.read(isize))

                print_matrix(
                        indices=indices,
                        mat_elems=elements,
                        no_elems=no_elem
                        )
                print()

                if no_elem < lembf:
                    end_of_file = True
                    break

            except:

                raise IOError('Error reading hamiltonian matrix elements...')
                exit()

    except:

        break

print("Total no. of non-zero elements read: ", elem_count)
