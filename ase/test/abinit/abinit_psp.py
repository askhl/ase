from ase.test import NotAvailable

try:
    from ase.test.abinit_installed import abinit_installed
    abinit_installed()
except NotAvailable:
    raise NotAvailable('Abinit required')

from ase.calculators.abinit import Abinit

from ase.utils.compound_test import MoleculeTest
molecule = MoleculeTest().compound

def get_run(atoms, pps, xc):
    ecut = 100
    kwargs = dict(
        label=xc + '_' + pps,
        xc=xc,
        width=0.001, # eV abinit does not accept 0.0
        ecut=ecut,
        pps=pps,
        )
    calc = Abinit(**kwargs)
    if pps == 'paw':
        calc.set_inp('pawecutdg', ecut) # pawecutdg should be ~ 3*ecut
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    return atoms


ref = {
    'PBE':
    {
    'fhi': -31.803974294357225,
    'hgh.k': -31.234649707555263,
    'hgh': -31.234277396315729, # MDTMP: uses LDA_HGH psp
    'tm': -32.15445642858279, #  MDTMP: uses LDA_TM psp
    'paw': -32.395413194274063,
    },
    'LDA':
    {
    'fhi': -31.265659502067816,
    'hgh.k': -30.712766446448665, # MDTMP: uses GGA_HGHK psp
    'hgh': -30.712407050439836,
    'tm': -31.593771991660827,
    'paw': -31.820373215999961,
    },
    }

formulas = ['H2']

pps = ['fhi', 'hgh', 'hgh.k', 'tm', 'paw']

energies = {}

for f in formulas:
    m = molecule(f)
    m.center(vacuum=1.5)
    for xc in ['LDA', 'PBE']:
        energies[xc] = {}
        for p in pps:
            try:
                atoms = get_run(m, p, xc)
                energies[xc][p] = atoms.get_potential_energy()
            except RuntimeError:
                pass # probably no pps found (or run failed)

# analyse results

for xc in ref.keys():
    for p in ref[xc].keys():
        diff = energies[xc][p] - ref[xc][p] # eV
        assert abs(diff) < 0.01, k + ': ' + str(diff)
