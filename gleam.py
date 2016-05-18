"""
Implementationg of the GLEaM model. See reference article:

Balcan, Duygu, Bruno Gonçalves, Hao Hu, José J. Ramasco, Vittoria Colizza,
and Alessandro Vespignani. 2010.
“Modeling the Spatial Spread of Infectious Diseases: The GLobal Epidemic and
Mobility Computational Model.”
Journal of Computational Science 1 (3): 132–45.
"""
# !/usr/bin/env python
#  encoding: utf-8


from copy import deepcopy
import math
import pandas as pd
import pygal
from pprint import pprint
import random
from utilities import plot_results
import numpy as np
import networkx as nx


class Subpopulation(dict):
    """
    Contains SLIR compartments for a given subpopulation.
    
    Attributes:
        id (TYPE): Description
        infectious_a (TYPE): Description
        infectious_nt (TYPE): Description
        infectious_t (TYPE): Description
        latent (TYPE): Description
        recovered (TYPE): Description
        susceptible (TYPE): Description
    """

    def __init__(self, _id, comps):
        """Summary
        
        Args:
            _id (TYPE): Description
            comps (TYPE): Description
        """
        self.id = _id
        self.susceptible = comps['susceptible']
        self.latent = comps['latent']
        self.infectious_nt = comps['infectious_nt']
        self.infectious_t = comps['infectious_t']
        self.infectious_a = comps['infectious_a']
        self.recovered = comps['recovered']

    def infectious_count(self):
        """Summary
        
        Returns:
            TYPE: Description
        """
        return self.infectious_nt + self.infectious_t + self.infectious_a

    def total_pop(self):
        """Summary
        
        Returns:
            TYPE: Description
        """
        return sum([self.susceptible, self.latent, self.infectious_nt,
                    self.infectious_t, self.infectious_a, self.recovered])
            

class Model(nx.DiGraph):
    """
    Uses nx.DiGraph to implement the network structure, and extends it with methods
    to run the simulations.
    
    Each node is a dict with values for each subpopulation compartment, ie
    susceptible, latent, infectious asymptomatic, infectious allowed to travel,
    infectious travel restricted, and recovered.
    
    Each edge is a commuting route between two nodes. Commuting rates between two
    subpopulations are encoded with the commuting_rate attribute for each edge.
    
    Attributes:
        commuting_return_rate (float): tau parameter - 1 / days spent out
        p_asymptomatic (float): probability of asymptomatic infection
        p_exit_latent (float): probability of exiting latent compartment
        p_recovery (float): probability of exiting infectious compartments
        p_travel_allowed (float): probability of being allowed to travel while infectious
    """
    def __init__(self, subpop_network, params):
        """
        DEFINE SUBPOP_NETWORK EXPECTED STRUCTURE
        
        Parameters
        ----------
        subpop_network : nx.DiGraph
            graph representation of the subpopulations and their commuting
            relationships
        params : dict
            model parameters
        
        Raises
        ------
        ValueError
        param dict does not contain all required values
        """
        super(Model, self).__init__(subpop_network)
        
        for param in ['p_exit_latent', 'p_recovery','p_asymptomatic',
                      'p_travel_allowed', 'commuting_return_rate']:
            if param not in params.keys():
                raise ValueError("Missing {} parameter".format(param))

        self.p_exit_latent = params['p_exit_latent']
        self.p_recovery = params['p_recovery']
        self.p_asymptomatic = params['p_asymptomatic']
        self.p_travel_allowed = params['p_travel_allowed']
        self.commuting_return_rate = params['commuting_return_rate']

    def coefficient_of_transmission(self, node_id):
        """Summary
        
        Args:
            node_id (int): id of the nx.DiGraph node
        
        Returns:
            float: Probability extracted from the binomial distribution
        """
        return np.random.binomial(self.node[node_id]['susceptible'],
                                  self.effective_force_of_infection(node_id))

    def transitions_from_latent(self, node_id):
        """
        First extract probabilities for asymptomatic, travelling and non travelling
        infectious among those that leave latent state, then scale by probability
        to exit latent state.
        
        Args:
            node_id (int): id of the nx.DiGraph node

        Returns:
            list of floats: Transition probabilies from latent to each infectious
                            compartments.
        """
        p_a = self.p_asymptomatic
        p_t = (1 - self.p_asymptomatic) * self.p_travel_allowed
        p_nt = (1 - self.p_asymptomatic) * (1 - self.p_travel_allowed)
        a, t, nt = np.random.multinomial(self.node[node_id]['latent'],
                                         [p_a, p_t, p_nt], size=1)
        return [x * self.p_exit_latent for x in [a, t, nt]]

    def effective_population(self, node_id):
        """
        For given subpopulation node, compute effective total population,
        taking into account the commuting rates between neighboring subpopulations.
        
        Args:
            node_id (int): id of the nx.DiGraph node

        Returns:
            int: Effective population of given node.
        """
        node_pop = self.node[node_id]
        node_exit_rate = self.get_exit_rate(node_id)
        local_pop = (node_pop['infectious_nt'] +
                     (sum(node_pop.values()) - node_pop.['infectious_nt']) /
                     (1 + node_exit_rate))

        other_pop = 0
        for nb_id in nx.neighbors(self, node_id):
            nb_pop = self.node[nb_id]
            nb_exit_rate = self.get_exit_rate(nb_id)
            other_pop += (((sum(nb_pop.values()) - nb_pop.infectious_nt) / 
                           (1 + nb_exit_rate)) *
                          self.edge[nb_id][node_id]['commuting_rate'] / self.return_rate)

        return local_pop + other_pop

    def effective_force_of_infection(self, node_id):
        """
        For given subpopulation node, compute effective for of infection,
        taking into account the commuting rates between neighboring subpopulations.
        
        Args:
            node_id (int): id of the nx.DiGraph node

        Returns:
            float: Effective force of infection at given node.
        """
        node_exit_rate = self.get_exit_rate(node_id)
        node_pop = self.node[node_id]
        
        foi_from_inside = 
        foi_from_outside = 

    def get_exit_rate(self, node_id):
        """
        Return the commuting exit rate for a given node, as a function of
        commuting and return rates.
        
        Args:
            node_id (int): id of the nx.DiGraph node
        """
        return (sum([e[2]['commuting_rate'] 
                    for e in self.out_edges(node_id, data=True)]) /
                self.return_rate)

    def seasonality(self, hemisphere):
        """
        Computes the scalar factor to apply on force of infection.
        
        Args:
            hemisphere (str): either 'north' or 'south'
        
        Returns:
            float: seasonality value for the given hemisphere
        """
        # TODO: implement seasonality computation
        return 1


    def total_infectious(self, node_pop):
        """Returns the total number of infectious people in a given node population.
        
        Args:
            node_pop (dict): Contains population values for the different compartments
        
        Returns:
            int: Sum of values in infectious compartments
        """
        return sum([node_pop.infectious_nt, node_pop.infectious_t, 
                    node_pop.infectious_a])
    def total_pop(self, node_pop):
        """Returns the total population in a given node population.
        
        Args:
            node_pop (dict): Contains population values for the different compartments

        Returns:
            int: Sum of values in all compartments
        """
        return sum(node_pop.values())


class CompartmentModel(object):
    """Summary
    
    Attributes:
        comps (TYPE): Description
        params (TYPE): Description
    """
    def __init__(self):
        """Summary
        """
        pass

    def compute_infection_rate(self):
        """Summary
        
        Returns:
            TYPE: Description
        """
        return (self.params['transmission_rate'] * self.infectious_count() /
                sum(self.comps.values()))

    def compute_commuting_rate

    def get_state(self):
        """
        Returns a copy of the model compartments.
        """
        return deepcopy(self.comps)

    def infectious_count(self):
        """Summary
        
        Returns:
            TYPE: Description
        """
        return (self.comps['symptomatic_travel'] +
                self.comps['symptomatic_no_travel'] +
                self.comps['asymptomatic'])

    def infect(self):
        """
        Updates counts in compartements according to model parameters.
        """

        inf_rate = self.compute_infection_rate()
        new_latent = math.ceil(self.comps['susceptible'] * inf_rate)
        new_sym_no_travel = math.ceil(self.comps['latent'] *
                                      self.params['p_exit_latent'] *
                                      (1 - self.params['p_asymptomatic_infection']) *
                                      (1 - self.params['p_travel_permission']))
        new_sym_travel = math.ceil(self.comps['latent'] *
                                   self.params['p_exit_latent'] *
                                   (1 - self.params['p_asymptomatic_infection']) *
                                   self.params['p_travel_permission'])
        new_asymptomatic = math.ceil(self.comps['latent'] *
                                     self.params['p_exit_latent'] *
                                     self.params['p_asymptomatic_infection'])
        new_snt_recovered = math.ceil(self.comps['symptomatic_no_travel'] *
                                      self.params['p_recovery'])
        new_st_recovered = math.ceil(self.comps['symptomatic_travel'] *
                                     self.params['p_recovery'])
        new_asym_recovered = math.ceil(self.comps['asymptomatic'] *
                                       self.params['p_recovery'])

        self.comps['susceptible'] -= new_latent
        self.comps['latent'] += new_latent
        self.comps['latent'] -= (new_sym_no_travel + new_sym_travel +
                                 new_asymptomatic)
        self.comps['symptomatic_no_travel'] += new_sym_no_travel
        self.comps['symptomatic_travel'] += new_sym_travel
        self.comps['asymptomatic'] += new_asymptomatic
        self.comps['symptomatic_no_travel'] -= new_snt_recovered
        self.comps['symptomatic_travel'] -= new_st_recovered
        self.comps['asymptomatic'] -= new_asym_recovered
        self.comps['recovered'] += (new_snt_recovered + new_st_recovered +
                                    new_asym_recovered)

    def set_populations(self, comps):
        """Summary
        
        Args:
            comps (TYPE): Description
        
        Returns:
            TYPE: Description
        
        Raises:
            Exception: Description
        """
        for x in ['susceptible', 'latent', 'symptomatic_no_travel',
                  'symptomatic_travel', 'asymptomatic', 'recovered']:
            if x not in comps.keys():
                raise Exception("GleamModel compartments must include {}".format(x))
        self.comps = comps

    def set_parameters(self, params):
        """Summary
        
        Args:
            params (TYPE): Description
        
        Returns:
            TYPE: Description
        
        Raises:
            Exception: Description
        """
        for x in ['p_exit_latent', 'transmission_rate', 'p_asymptomatic_infection',
                  'p_travel_permission', 'p_recovery']:
            if x not in params.keys():
                raise Exception("GleamModel parameters must include {}".format(x))
        self.params = params

if __name__ == '__main__':
    parameters = {'p_exit_latent': 0.1,
                  'transmission_rate': 0.5,
                  'p_asymptomatic_infection': 0.1,
                  'p_travel_permission': 0.3,
                  'p_recovery': 0.05}

    populations = {'susceptible': 1950,
                    'latent': 0,
                    'symptomatic_no_travel': 10,
                    'symptomatic_travel': 0,
                    'asymptomatic': 0,
                    'recovered': 0}
    gm = CompartmentModel()
    gm.set_populations(populations)
    gm.set_parameters(parameters)
    results = []
    for i in range(1000):
        results.append(gm.get_state())
        gm.infect()
    save_file = '/home/hugo/Projects/networks/network_theory_application_course_project/visualization/svg_test.svg'
    print
    plot_results(results)
