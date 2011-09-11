def get_bandgap(atoms, mode='direct', spin=0):
    import numpy as np
    assert mode in ['direct', 'indirect']
    calculator = atoms.get_calculator()
    direct_bandgap = None
    indirect_bandgap = None
    nelect = calculator.get_number_of_electrons()
    nband = calculator.get_number_of_bands()
    magmom = calculator.get_magnetic_moment(atoms)
    occupations = np.array(calculator.get_occupation_numbers(spin=spin))
    nhomo = len([a for a in occupations if a > 0.0]) - 1
    nlumo = nhomo + 1
    assert nband >= nlumo + 1
    # find direct and indirect (makes no sense for molecules) band gap
    ehomo = [calculator.get_eigenvalues(kpt=0, spin=spin)[nhomo]]
    elumo = [calculator.get_eigenvalues(kpt=0, spin=spin)[nlumo]]
    direct_bandgap = elumo[0] - ehomo[0]
    ibz_k_points = calculator.get_ibz_k_points()
    for k in range(1, len(ibz_k_points)):
        ehomo.append(calculator.get_eigenvalues(kpt=k, spin=spin)[nhomo])
        elumo.append(calculator.get_eigenvalues(kpt=k, spin=spin)[nlumo])
        direct_bandgap = min(direct_bandgap, elumo[k]-ehomo[k])
    if len(ibz_k_points) > 1:
        indirect_bandgap =  min([elumo[k] for k in range(len(elumo))]) - max([ehomo[k] for k in range(len(ehomo))])
    if mode == 'direct':
        return direct_bandgap
    elif mode == 'indirect':
        return indirect_bandgap
