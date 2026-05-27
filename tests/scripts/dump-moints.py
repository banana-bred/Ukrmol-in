#!/usr/bin/env python3
#
# Convert `moints` to text
# ~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 dump-moints.py
#
# Writes subset of the content of the `moints` file (molecular integrals) to standard output
# for debugging. No arguments are needed. The file `moints` is expected to be
# in the working directory.
#
# Assumes that the file uses int64 and real64 types.
#

from array  import array
from struct import unpack

# The types below need to be compatible with GBTOlib INTEGERs and REALs.
isize = 4; iform = 'i'
lsize = 8; lform = 'l'
rsize = 8; rform = 'd'
#rsize = 16; rform = 'q'

line_len = 207
file_format = ''
format_10 = 'DATA FILE version 1.0'
format_11 = 'DATA FILE version 1.1'

# Print array of real values as table
def tabulate (frm, ncols, data):
    nrows = len(data) // ncols + 1
    for i in range(0, nrows):
        a = i * ncols
        b = min((i + 1)*ncols, len(data))
        for j in range(a, b):
            print(frm.format(data[j]), end = '')
        print()


# Read array of the given number of real numbers from file
def read_reals(f, n):
    if (rform == 'q'):
        f.read(n*rsize)  # advance file pointer but ignore the data
        return n*('(quad)',)
    else:
        return unpack(rform*n, f.read(n*rsize))


def read_ijkl_indices(f):
    nints, ntypes = unpack(lform * 2, f.read(2 * lsize))
    print('- nints = {}'.format(nints))
    print('- ntypes = {}'.format(ntypes))
    ints = array(lform); ints.fromfile(f, nints*ntypes)
    print('- integral indices')
    tabulate('{:5d}', 30, ints)


def read_1el_integrals(f):
    prec, tol, a = read_reals(f, 3)
    print('- prec = {}'.format(prec))
    print('- tol = {}'.format(tol))
    print('- Ra = {}'.format(a))
    max_ijrs_size, = read_reals(f, 1)
    two_p_continuum, use_spherical_cgto_alg, max_property_l, calculate_overlap_ints, calculate_kinetic_energy_ints, \
    calculate_nuclear_attraction_ints, calculate_one_el_hamiltonian_ints, calculate_property_ints, calculate_two_el_ints \
        = unpack(lform * 9, f.read(9 * lsize))
    dipole_damp_factor, = read_reals(f, 1) if file_format >= format_11 else (0.0,)
    print('- max_ijrs_size = {}'.format(max_ijrs_size))
    print('- two_p_continuum = {}'.format(two_p_continuum))
    print('- use_spherical_cgto_alg = {}'.format(use_spherical_cgto_alg))
    print('- max_property_l = {}'.format(max_property_l))
    print('- calculate_overlap_ints = {}'.format(calculate_overlap_ints))
    print('- calculate_kinetic_energy_ints = {}'.format(calculate_kinetic_energy_ints))
    print('- calculate_nuclear_attraction_ints = {}'.format(calculate_nuclear_attraction_ints))
    print('- calculate_property_ints = {}'.format(calculate_property_ints))
    print('- calculate_two_el_ints = {}'.format(calculate_two_el_ints))
    print('- dipole_damp_factor = {}'.format(dipole_damp_factor))
    m, n = unpack(lform * 2, f.read(2 * lsize))
    print('- array dimensions: {} x {}'.format(m, n))
    for i in range(0, n):
        column_descriptor = bytes(unpack(line_len * 'b', f.read(line_len))).decode('utf-8').strip()
        print('- column #{} descriptor: "{}"'.format(i + 1, column_descriptor))
    for i in range(0, n):
        vals = array('d')
        vals.fromfile(f, m*rsize//8)
        if (rform == 'd'):
            print('- integral values for column #{}'.format(i + 1))
            tabulate('{:25.15e}', 6, vals)


def read_2el_integrals(f):
    prec, tol, a = read_reals(f, 3)
    print('- prec = {}'.format(prec))
    print('- tol = {}'.format(tol))
    print('- Ra = {}'.format(a))
    max_ijrs_size, = read_reals(f, 1)
    two_p_continuum, use_spherical_cgto_alg, max_property_l, calculate_overlap_ints, calculate_kinetic_energy_ints, \
    calculate_nuclear_attraction_ints, calculate_one_el_hamiltonian_ints, calculate_property_ints, calculate_two_el_ints \
        = unpack(lform * 9, f.read(9 * lsize))
    dipole_damp_factor, = read_reals(f, 1) if file_format >= format_11 else (0.0,)
    print('- max_ijrs_size = {}'.format(max_ijrs_size))
    print('- two_p_continuum = {}'.format(two_p_continuum))
    print('- use_spherical_cgto_alg = {}'.format(use_spherical_cgto_alg))
    print('- max_property_l = {}'.format(max_property_l))
    print('- calculate_overlap_ints = {}'.format(calculate_overlap_ints))
    print('- calculate_kinetic_energy_ints = {}'.format(calculate_kinetic_energy_ints))
    print('- calculate_nuclear_attraction_ints = {}'.format(calculate_nuclear_attraction_ints))
    print('- calculate_property_ints = {}'.format(calculate_property_ints))
    print('- calculate_two_el_ints = {}'.format(calculate_two_el_ints))
    print('- dipole_damp_factor = {}'.format(dipole_damp_factor))
    m, n = unpack(lform * 2, f.read(2 * lsize))
    print('- array dimensions: {} x {}'.format(m, n))
    for i in range(0, n):
        column_descriptor = bytes(unpack(line_len * 'b', f.read(line_len))).decode('utf-8').strip()
        print('- column #{} descriptor: "{}"'.format(i + 1, column_descriptor))
    for i in range(0, n):
        vals = array('d')
        vals.fromfile(f, m*rsize//8)
        if (rform == 'd'):
            print('- integral values for column #{}'.format(i + 1))
            tabulate('{:25.15e}', 6, vals)


# Open the file for binary reading.
f = open('moints', 'rb')

print('File "moints"')

file_format = bytes(unpack(line_len * 'b', f.read(line_len))).decode('utf-8').strip()
print('- title: "{}"'.format(file_format))

flt_precision, int_bit_size = unpack(iform * 2, f.read(2 * isize))
print('- flt_precision = {}'.format(flt_precision))
print('- int_bit_size = {}'.format(int_bit_size))

no_headers, = unpack(lform, f.read(lsize))
print('- no_headers = {}'.format(no_headers))

for ih in range(0, no_headers):
    
    print('\nHeader #{}'.format(ih + 1))

    header_title = bytes(unpack(line_len * 'b', f.read(line_len))).decode('utf-8').strip()
    print('- title: "{}"'.format(header_title))

    first_record, last_record, next_data = unpack(lform * 3, f.read(3 * lsize))
    print('- first_record = {}'.format(first_record))
    print('- last_record = {}'.format(last_record))
    print('- next_data = {}'.format(next_data))
    
    if header_title == 'molecular_orbital_basis_obj:One-particle AO->MO transformed integrals for symmetric operators':
        read_1el_integrals(f)
    
    if header_title == 'molecular_orbital_basis_obj:Two-particle AO->MO transformed integrals for symmetric operators':
        read_2el_integrals(f)
        
    if header_title == 'molecular_orbital_basis_obj:ijkl indices':
        read_ijkl_indices(f)
    
    if next_data > 0:
        f.seek(next_data - 1)
