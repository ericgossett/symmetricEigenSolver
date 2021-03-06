"""
Our example run of pyQuante 1.6.5

Note the resulting orbital energies and molecular orbitals are those of the
final iteration of the SCF calculation (see ham.out files for the output
of each iteration)
"""
from PyQuante.Molecule import Molecule
from PyQuante.dft import dft
from PyQuante.hartree_fock import rhf

h2 = Molecule(
  'H2',
  [
    (1,  (0.00000000,     0.00000000,     0.36628549)),
    (1,  (0.00000000,     0.00000000,    -0.36628549))
   ],
  units='Angstrom'
)

h2o = Molecule(
    'H2O',
    [
        (8, (0.00000000,     0.00000000,     0.04851804)),
        (1, (0.75300223,     0.00000000,    -0.51923377)),
        (1, (-0.75300223,     0.00000000,    -0.51923377))
    ],
    units='Angstrom'
)

benz = Molecule(
    'benzene',
    [
        (6,  (0.000,  1.396,  0.000)),
        (6,  (1.209,  0.698,  0.000)),
        (6,  (1.209, -0.698,  0.000)),
        (6,  (0.000, -1.396,  0.000)),
        (6, (-1.209, -0.698,  0.000)),
        (6, (-1.209,  0.698,  0.000)),
        (1,  (0.000,  2.479,  0.000)),
        (1,  (2.147,  1.240,  0.000)),
        (1,  (2.147, -1.240,  0.000)),
        (1,  (0.000, -2.479,  0.000)),
        (1, (-2.147, -1.240,  0.000)),
        (1, (-2.147,  1.240,  0.000))
    ],
    units='Angstrom'
)

print "------------     h2      ------------"
en, orbe, orbs = rhf(h2)

print 'orbe', orbe
print 'orbs', orbs

print "------------     h2o     ------------"

en, orbe, orbs = rhf(h2o)

print 'orbe', orbe
print 'orbs', orbs

print "------------   beneze   ------------"

en, orbe, orbs = rhf(benz)

print 'orbe', orbe
print 'orbs', orbs
