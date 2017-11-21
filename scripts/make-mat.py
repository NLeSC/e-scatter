#!/usr/bin/env python3
import os
import glob
import subprocess
import sys
import numpy as np
import scipy.integrate as integrate
from constants import (N_A, eV, h, pi, k_B, h_bar, T)

work_dir = os.getcwd()
elsepa_dir = os.path.join(work_dir, "data/elsepa")
endf_dir = os.path.join(work_dir, "data/ENDF-B-VII.1")

def run_elsepa(Z, out_fn):
    out_xml = open(out_fn, 'w')
    os.chdir(elsepa_dir)
    for f in glob.glob('*.dat'):
        os.remove(f)

    no_muffin_Z = [1, 7, 8]

    mexch = 1
    mcpol = 2
    mabs = 1
    muffin = 1
    if Z in no_muffin_Z:
        muffin = 0

    elscata_in = ''
    # atomic number                             [none]
    elscata_in += 'IZ      {}\n'.format(Z)
    # rho_n (1=P, 2=U, 3=F, 4=Uu)                [  3]
    elscata_in += 'MNUCL   {}\n'.format(3)
    # rho_e (1=TFM, 2=TFD, 3=DHFS, 4=DF, 5=file) [  4]
    elscata_in += 'MELEC   {}\n'.format(4)
    # 0=free atom, 1=muffin-tin model            [  0]
    elscata_in += 'MUFFIN  {}\n'.format(muffin)
    # -1=electron, +1=positron                   [ -1]
    elscata_in += 'IELEC   {}\n'.format(-1)
    # V_ex (0=none, 1=FM, 2=TF, 3=RT)            [  1]
    elscata_in += 'MEXCH   {}\n'.format(mexch)
    # V_cp (0=none, 1=B, 2=LDA)                  [  0]
    elscata_in += 'MCPOL   {}\n'.format(mcpol)
    # high-E factorization (0=no, 1=yes, 2=Born) [  1]
    elscata_in += 'IHEF    {}\n'.format(0)
    # W_abs (0=none, 1=LDA)                      [  0]
    elscata_in += 'MABS    {}\n'.format(mabs)

    for E in [10, 20, 30, 40, 50, 60, 70, 80, 90,
              100, 200, 300, 400, 500, 600, 700, 800, 900,
              1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
              10000]:
        # kinetic energy (eV)                       [none]
        elscata_in += 'EV      {}\n'.format(E)

    subprocess.run(
        [os.path.join(elsepa_dir, 'elscata')],
        input=bytes(elscata_in, 'UTF-8'), stdout=subprocess.DEVNULL, check=True)

    os.chdir(work_dir)

    subprocess.run(
        ["scripts/elsepa-parse.py", elsepa_dir],
        stdout=out_xml, check=True)


def run_phonon(eps_ac, c_s, rho_m, M_tot, lattice, out_fn):
    out_xml = open(out_fn, 'w')
    subprocess.run(
        ['scripts/phonon-gen.py', '--eps_ac', str(eps_ac), '--c_s', str(c_s),
         '--rho_m', str(rho_m), '--M', str(M_tot), '--a', str(lattice)],
        stdout=out_xml, check=True)


def run_endf(Z, out_fn):
    out_xml = open(out_fn, 'w')
    subprocess.run(
        ['scripts/endf-parse.py', endf_dir, str(Z)],
        stdout=out_xml, check=True)

def phonon_loss(c_s, a, T):
    w = lambda k: 2*c_s/a*np.sin(k*a/2)
    N = lambda k: 1/np.expm1(h_bar*w(k)/(k_B*T))

    x1 = lambda k: w(k)*k*k
    x2 = lambda k: (2*N(k)+1)*k*k
    y1, err1 = integrate.quad(x1, 0, pi/a)
    y2, err2 = integrate.quad(x2, 0, pi/a)

    E_ph = h_bar * y1/y2
    return E_ph

if len(sys.argv) != 2:
    print("usage: {} <material_definition.py>".format(sys.argv[0]))
    sys.exit(1)

mat_def_fn = sys.argv[1]
with open(mat_def_fn) as f:
    mat_def = f.read()
mat = eval(mat_def)

mat['M_tot'] = 0
for elem in mat['elements']:
    mat['M_tot'] += mat['elements'][elem]['count']*mat['elements'][elem]['M']
mat['rho_n'] = N_A/mat['M_tot'] * mat['rho_m']
mat['phonon_loss'] = phonon_loss(mat['c_s'], mat['lattice'], T)

mat_dir = os.path.join(work_dir, "data/materials", mat['name'])
el_phon_fn = os.path.join(mat_dir, 'elastic-phonon.xml')
el_mott_fn = os.path.join(mat_dir, 'elastic-mott.xml')
el_fn = os.path.join(mat_dir, 'elastic.xml')
inel_fn = os.path.join(mat_dir, 'inelastic.xml')
ion_fn = os.path.join(mat_dir, 'ionization.xml')
osi_fn = os.path.join(mat_dir, 'outer_shell.dat')

with open(el_mott_fn, 'w') as el_mott_xml:
    el_mott_xml.write('<cstable type="elastic">\n</cstable>\n')

for elem in mat['elements']:
    Z = mat['elements'][elem]['Z']
    count = mat['elements'][elem]['count']

    element_el_mott_fn = os.path.join(
        mat_dir, 'elastic-mott-{}.xml'.format(elem))
    print('running Elsepa for element {}'.format(elem), file=sys.stderr)
    run_elsepa(Z, element_el_mott_fn)
    with open('tmp.xml', 'w') as tmp_xml:
        subprocess.run(
            ['bin/cstool', 'mad', str(1), el_mott_fn, str(count),
             element_el_mott_fn],
             stdout=tmp_xml, check=True)
    os.rename('tmp.xml', el_mott_fn)

    element_ion_fn = os.path.join(mat_dir, 'ionization-{}.xml'.format(elem))
    print('parsing ENDF files for element {}'.format(elem), file=sys.stderr)
    run_endf(Z, element_ion_fn)
    with open(ion_fn, 'w') as ion_xml:
        subprocess.run(
            ['bin/cstool', 'mad', str(count), element_ion_fn],
            stdout=ion_xml, check=True)

print('generating phonon cross-sections', file=sys.stderr)
run_phonon(mat['eps_ac'], mat['c_s'], mat['rho_m'],
    mat['M_tot'], mat['lattice'], el_phon_fn)

with open(el_fn, 'w') as el_xml:
    subprocess.run(
        ['bin/cstool', 'merge', el_phon_fn, el_mott_fn, '100*eV', '200*eV'],
        stdout=el_xml, check=True)

with open(mat['elf-file'], 'r') as elf_f, open(osi_fn, 'w') as osi_f:
    osi_f.write(elf_f.readline())

print('generating inelastic cross-sections', file=sys.stderr)
with open(inel_fn, 'w') as inel_xml:
    subprocess.run(
        ['scripts/inelastic-gen.py', '--number-density', str(mat['rho_n']),
         '--elf-file', mat['elf-file'], '--fermi', str(mat['fermi'] / eV)],
        stdout=inel_xml, check=True)

print('compiling material', file=sys.stderr)
process = ['bin/cstool', 'compile-mat',
     '{}/{}.mat'.format(mat_dir, mat['name']),
     '--name', mat['name'],
     '--elastic', el_fn,
     '--inelastic', inel_fn,
     '--ionization', ion_fn,
     '--outer-shell', osi_fn,
     '--number-density', str(mat['rho_n']),
     '--fermi-energy', str(mat['fermi']),
     '--work-function', str(mat['work_func']),
     '--phonon-loss', str(mat['phonon_loss'])]
if 'band_gap' in mat:
    process.extend(('--band-gap', str(mat['band_gap'])))
subprocess.run(process, check=True)
