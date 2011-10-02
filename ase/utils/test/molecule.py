from ase.utils.test.compound import CompoundTestEnergy

class MoleculeTestEnergy(CompoundTestEnergy):

    def __init__(self, calculator, calculate=None, retrieve=None,
                 name='test', data=None,
                 vacuum=8.0):

        CompoundTestEnergy.__init__(self, calculator, calculate, retrieve,
                                    name=name, data=data)

        self.vacuum = vacuum

    def setup_system(self, formula):
        """Create an Atoms object for the given formula.

        By default this will be loaded from the data, setting
        the cell size by means of the molecule test's vacuum parameter."""
        system = self.data[formula]
        if self.vacuum is not None:
            system.center(vacuum=self.vacuum)
        return system

def get_atomization_reactions(data):
    """Get definitions of atomization reactions for the database data.  """
    reactions = []
    formulas = data.keys()
    formulas.sort()
    for formula in formulas:
        if len(data[formula]) == 1:
            continue # skip atom
        else:
            stoich = []
            # find constituing atoms in the database
            for a in data[formula].get_chemical_symbols():
                for f in formulas:
                    fsymbols = data[f].get_chemical_symbols()
                    if len(fsymbols) == 1 and fsymbols[0] == a:
                        stoich.append((f, 1))
            stoich.append((formula, - 1))
            # for atomization reactions use formula name as reaction identifier
            stoich.append(('reaction_id', formula))
            reactions.append(stoich)
    return reactions
