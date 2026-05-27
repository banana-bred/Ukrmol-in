#!/usr/bin/env python3
#
# Visualize CONGEN output as plot (optionally save all to pngs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Typical usage:
#
#     python3 congen_visualizer.py fort.70
#
#     (it will ask user to select the CSF they wish to print)
#     (If user enters '0' then all pngs will be saved)
#
# Optional usage: user can supply 'core' orbitals which helps to differentiate
#                 between core and virtual orbitals but it is not necessary
#                 e.g. for c2v symmtery with core = [5,1,1,0]
#
#     python3 congen_visualizer.py fort.70 5,1,1,0
#
# Assumes that the file uses int64 and real64 types.
#
# It can also be useful to import this script into a jupyter notebook
# e.g. from congen_visualizer import read_congen
# then you can see the pictoral representation of several determinants
# all in one place.
#

from array  import array
from struct import unpack
from sys    import argv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

linewidth = 1

class CSF:

    def __init__(self, dets, coeffs, id):
        self.determinants = dets
        self.coefficients = coeffs
        self.id = id
        self.no_dets = len(dets)

    def __str__(self):

        string = f'\n--- displaying CSF {self.id}: ---'
        string += f'\nno of dets: {self.no_dets}'
        string += f'\ncoefficients : {self.coefficients}\n'
        for i, det in enumerate(self.determinants):
            string += f'\n--- determinant #{i+1}: ---\n'
            string += str(det)
            string += '\n'

        return string

    def visualize(self, figsize=(3, 3)):

        fig, axes = plt.subplots(
                1,
                self.no_dets,
                figsize=(figsize[0]*self.no_dets, figsize[1]),
                )

        if self.no_dets == 1:
            axes = [axes]

        for i, det in enumerate(self.determinants):
            det.visualize(ax=axes[i])
            title = f'CSF#{self.id:3d} c={self.coefficients[i]:5.4f}'
            axes[i].set_title(title)

        return fig

class Determinant:

    def __init__(self, reference, orb_table):
        self.reference = reference
        self.determinant = reference
        self.substitutions = dict()
        self.orb_table = orb_table

    def substitute(self, orb_a, orb_b):
        substitutions = dict(zip(orb_a, orb_b))
        self.substitutions = substitutions

        temp_det = []

        for orb in self.reference:
            try:
                temp_det.append(substitutions[orb])
            except KeyError:
                temp_det.append(orb)


        self.determinant = temp_det
        pass

    def __str__(self):
        string = f'reference    : {self.reference}'
        string += f'\ndeterminant  : {self.determinant}'
        string += f'\norbital subs :'
        for i, f in self.substitutions.items():
            string += f'\n               ({i} -> {f})'
        return string

    def visualize(self, ax=None):

        if ax is None:
            figsize = (3, 3)
            ax = plt.figure(figsize=figsize).gca()

        df = self.orb_table
        max_core = df[df['group']=='core']['orbital'].max() + 1
        max_virt = df[df['group']=='virt']['orbital'].max() + 1
        max_sym = df['symmetry'].max() + 1

        if np.isnan(max_core):
            max_core = 0


        draw_guides(ax=ax, core=max_core, virt=max_core+max_virt, max_x=4)

        for orb in self.determinant:
            orbital, symmetry, spin, group = df.loc[orb].values
            if group == 'core':
                draw_orbital(ax=ax, num=orbital, sym=symmetry, offset=0)
                draw_electron(ax=ax, num=orbital, sym=symmetry, spin=spin, offset=0)
            elif group == 'virt':
                draw_orbital(ax=ax, num=orbital, sym=symmetry, offset=max_core)
                draw_electron(ax=ax, num=orbital, sym=symmetry, spin=spin, offset=max_core)
            elif group == 'cont':
                draw_orbital(ax=ax, num=orbital, sym=symmetry, offset=max_core+max_virt)
                draw_electron(ax=ax, num=orbital, sym=symmetry, spin=spin, offset=max_core+max_virt)
            else:
                print(f'ERROR :: group {group} not recognized.')

        ax.set_xlim([0, max_sym])
        ax.set_ylim([-0.5, max_core+max_virt+3])
        ax.set_xticks(np.arange(0.5, max_sym, 1))
        ax.set_xticklabels(list(range(max_sym)))
        ax.get_yaxis().set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)

        pass

def draw_guides(ax, core, virt, max_x):
    if core > 0:
        line = plt.Line2D((0, max_x), (core, core), lw=linewidth, ls='--', color='b')
        ax.add_line(line)
    line = plt.Line2D((0, max_x), (virt, virt), lw=linewidth, ls='--', color='r')
    ax.add_line(line)
    pass

def draw_orbital(ax, num, sym, offset):
    xbuf = 0.1
    mo_width = 1.0
    line = plt.Line2D((sym+xbuf, sym+mo_width-xbuf), (num+offset, num+offset), lw=linewidth, ls='-', color='k')
    ax.add_line(line)
    pass

def draw_electron(ax, num, sym, spin, offset):
    electron_height = 1.0
    ybuf = 0.1

    if spin == 0: # spin up

        ax.arrow(sym+0.4, num+ybuf+offset, dx=0, dy=electron_height-2*ybuf,
                       shape='right',
                       length_includes_head=True,
                       head_width=0.1,
                       head_length=0.1,
                       lw=linewidth, 
                       ls='-',
                       color='k',
                      )

    else: # spin down

        ax.arrow(sym+0.6, num+electron_height-ybuf+offset, dx=0, dy=-electron_height+2*ybuf,
                       shape='right',
                       length_includes_head=True,
                       head_width=0.1,
                       head_length=0.1,
                       lw=linewidth, 
                       ls='-',
                       color='k',
                      )

    pass

def generate_orbital_table(nob0, nob, core=None):

    #nob0 = target orbs per irrep.
    #nob  = no. of orbs per irrep.
    #core = no. of core orbs per irrep. (optional - can help with visualization)

    nsym = len(nob)
    if core is None:
        core = [0]*nsym
    group = ''
    spin_orbital_table = []
    orb_id = 0

    for sym in range(nsym):
        for orb in range(nob[sym]):
            for spin in [0, 1]:
                if orb < core[sym]:
                    group = 'core'
                    orb_id = orb
                elif orb < nob0[sym]:
                    group = 'virt'
                    orb_id = orb - core[sym]
                else:
                    group = 'cont'
                    orb_id = orb - nob0[sym]
                spin_orbital_table.append([orb_id, sym, spin, group])

    spin_orbital_table = pd.DataFrame(spin_orbital_table)
    spin_orbital_table.columns = ['orbital', 'symmetry', 'spin', 'group']

    spin_orbital_table.index = spin_orbital_table.index + 1

    return spin_orbital_table

def read_congen(filename, core=None):

    csfs = []

    # Expect electron (= 0) or positron (= 1) configuration set.
    iposit = 0

    # spin orbital table for determinants
    spin_orbital_table = None

    # The types below need to be compatible with CONGEN INTEGERs and REALs.
    isize = 4; iform = 'i'
    lsize = 8; lform = 'l'
    rsize = 8; rform = 'd'

    # Open the file for binary reading.
    f = open(filename, 'rb')

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
            else:
                begin, = unpack(iform, f.read(isize))
                iphz = array(lform); iphz.fromfile(f, nctarg)
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

            spin_orbital_table = generate_orbital_table(nob0=nob0, nob=nob, core=core)

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
        i = indo[icsf]-1
        dets = []
        coeffs = []
        for idet in range(icdo[icsf], icdo[icsf+1]):
            nr = ndo[i]
            src = [j for j in ndo[i+1:i+nr+1]]
            dst = [j for j in ndo[i+nr+1:i+2*nr+1]]
            i = i + 2*nr + 1
            det = Determinant(reference=ndtrf.tolist(), orb_table=spin_orbital_table)
            det.substitute(orb_a=src , orb_b=dst)
            dets.append(det)
            coeffs.append(cdo[idet-1])

        csfs.append(CSF(dets=dets, coeffs=coeffs, id=icsf + 1))

    return csfs


def main():

    def repeatOnError(*exceptions):
        def checking(function):
            def checked(*args, **kwargs):
                while True:
                    try:
                        result = function(*args, **kwargs)
                        if result > ncsf:
                            raise ValueError
                        if result < 0:
                            raise ValueError
                    except exceptions as problem:
                        print("There was a problem with the input.")
                        print("Please enter valid integer.")
                    else:
                        return result
            return checked
        return checking

    @repeatOnError(ValueError)
    def get_CSF_id(ncsf):
        msg = f'Please enter CSF number (1-{ncsf} or press 0 to save all CSFs as PNGs): '
        return int(input(msg))

    filename = argv[1]
    try:
        core = argv[2]
        core = [int(i) for i in core.split(',')]
    except:
        print('optional core orbitals not supplied.')
        core = None

    csfs = read_congen(filename=filename, core=core)
    ncsf = len(csfs)

    csf_id = get_CSF_id(ncsf)

    if csf_id > 0:
        csf = csfs[csf_id-1]
        fig = csf.visualize()
        print(csf)
        plt.show()
        fig.savefig(f'csf{csf_id:05d}.png')

    else:
        for i, csf in enumerate(csfs):
            fig = csf.visualize()
            fig.savefig(f'csf{i+1:05d}.png')
            plt.close()

    return 0

if __name__ == "__main__":
    main()
