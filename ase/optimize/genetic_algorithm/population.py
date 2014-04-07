""" Implementaiton of a population for maintaining a GA population and
proposing structures to pair. """
from random import randrange, random
from math import tanh, sqrt


def count_looks_like(a, all_cand, comp):
    """ Utility method for counting occurences. """
    n = 0
    for b in all_cand:
        if a.info['confid'] == b.info['confid']:
            continue
        if comp.looks_like(a, b):
            n += 1
    return n


class Population(object):
    """
       Population class which maintains the current population
       and proposes which candidates to pair together.

       Parameters:

       data_connection: DataConnection object
       population_size: The number of candidates in the population
       comparator: Comparator object which can tell if two
       configurations are equal.
    """
    def __init__(self, data_connection, population_size, comparator):
        self.dc = data_connection
        self.pop_size = population_size
        self.comparator = comparator
        self.pop = []
        self.pairs = None
        self.all_cand = None
        self.__initialize_pop__()

    def __initialize_pop__(self):
        """ Private method that initalizes the population when
            the population is created. """

        # Get all relaxed candidates from the database
        all_cand = self.dc.get_all_relaxed_candidates()
        all_cand.sort(key=lambda x: x.get_potential_energy())

        # Fill up the population with the self.pop_size most stable
        # unique candidates.
        i = 0
        while i < len(all_cand) and len(self.pop) < self.pop_size:
            c = all_cand[i]
            i += 1
            eq = False
            for a in self.pop:
                if self.comparator.looks_like(a, c):
                    eq = True
                    break
            if not eq:
                self.pop.append(c)

        for a in self.pop:
            a.info['looks_like'] = count_looks_like(a, all_cand,
                                                    self.comparator)

        self.all_cand = all_cand
        self.__calc_participation__()

    def __calc_participation__(self):
        """ Determines, from the database, how many times each
            candidate has been used to generate new candidates. """
        (participation, pairs) = self.dc.get_participation_in_pairing()
        for a in self.pop:
            if a.info['confid'] in participation.keys():
                a.info['n_paired'] = participation[a.info['confid']]
            else:
                a.info['n_paired'] = 0
        self.pairs = pairs

    def update(self):
        """ New candidates can be added to the database
            after the population object has been created.
            This method extracts these new candidates from the
            database and includes them in the population. """

        if len(self.pop) == 0:
            self.__initialize_pop__()

        new_cand = self.dc.get_all_relaxed_candidates(only_new=True)
        for a in new_cand:
            self.__add_candidate__(a)
            self.all_cand.append(a)
        self.__calc_participation__()

    def get_current_population(self):
        """ Returns a copy of the current population. """
        self.update()
        return [a.copy() for a in self.pop]

    def __add_candidate__(self, a):
        """ Adds a single candidate to the population. """

        #check if the structure is too high in energy
        if a.get_potential_energy() > self.pop[-1].get_potential_energy() \
                and len(self.pop) == self.pop_size:
            return

        # check if the new candidate should
        # replace a similar structure in the population
        for (i, b) in enumerate(self.pop):
            if self.comparator.looks_like(a, b):
                if b.get_potential_energy() > a.get_potential_energy():
                    del self.pop[i]
                    a.info['looks_like'] = count_looks_like(a,
                                                            self.all_cand,
                                                            self.comparator)
                    self.pop.append(a)
                    self.pop.sort(key=lambda x: x.get_potential_energy())
                return

        # the new candidate needs to be added, so remove the highest
        # energy one
        if len(self.pop) == self.pop_size:
            del self.pop[-1]

        # add the new candidate
        a.info['looks_like'] = count_looks_like(a,
                                                self.all_cand,
                                                self.comparator)
        self.pop.append(a)
        self.pop.sort(key=lambda x: x.get_potential_energy())

    def __get_fitness__(self, indecies, with_history=True):
        """ Calculates the fitness using the formula from
             L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
        """

        energies = [x.get_potential_energy() for x in self.pop]
        min_e = min(energies)
        max_e = max(energies)
        T = max_e - min_e
        if isinstance(indecies, int):
            indecies = [indecies]

        f = [0.5 * (1. - tanh(2. * (energies[i] - min_e) / T - 1.))
             for i in indecies]
        if with_history:
            M = [float(self.pop[i].info['n_paired']) for i in indecies]
            L = [float(self.pop[i].info['looks_like']) for i in indecies]
            f = [f[i] * 1. / sqrt(1. + M[i]) * 1. / sqrt(1. + L[i])
                 for i in xrange(len(f))]
        return f

    def get_two_candidates(self, with_history=True):
        """ Returns two candidates for pairing employing the
            fitness criteria from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            and the roulete wheel selection scheme described in
            R.L. Johnston Dalton Transactions,
            Vol. 22, No. 22. (2003), pp. 4193-4207
        """

        if len(self.pop) < 2:
            self.update()

        if len(self.pop) < 2:
            return None

        fit = self.__get_fitness__(range(len(self.pop)), with_history)
        fmax = max(fit)
        c1 = self.pop[0]
        c2 = self.pop[0]
        used_before = False
        while c1.info['confid'] == c2.info['confid'] and not used_before:
            nnf = True
            while nnf:
                t = randrange(0, len(self.pop), 1)
                if fit[t] > random() * fmax:
                    c1 = self.pop[t]
                    nnf = False
            nnf = True
            while nnf:
                t = randrange(0, len(self.pop), 1)
                if fit[t] > random() * fmax:
                    c2 = self.pop[t]
                    nnf = False

            c1id = c1.info['confid']
            c2id = c2.info['confid']
            used_before = (min([c1id, c2id]), max([c1id, c2id])) in self.pairs
        return (c1.copy(), c2.copy())
