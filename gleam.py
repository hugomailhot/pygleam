# !/usr/bin/env python
#  encoding: utf-8

"""
Implementationg of the GLEaM model. See reference article:

Balcan, Duygu, Bruno Gonçalves, Hao Hu, José J. Ramasco, Vittoria Colizza,
and Alessandro Vespignani. 2010.
“Modeling the Spatial Spread of Infectious Diseases: The GLobal Epidemic and
Mobility Computational Model.”
Journal of Computational Science 1 (3): 132–45.
"""

from copy import deepcopy
import math
import pandas as pd
import pygal
from pprint import pprint
import random
from utilities import plot_results
import networkx as nx


class Subpopulation(object):
    """
    Contains SLIR compartments for a given subpopulation.
    """

    def __init__(self, _id, comps):
        self.id = _id
        self.susceptible = comps['susceptible']
        self.latent = comps['latent']
        self.infectious_nt = comps['infectious_nt']
        self.infectious_t = comps['infectious_t']
        self.infectious_a = comps['infectious_a']
        self.recovered = comps['recovered']

        self.effective_force_of_infection = 0
        self.effective_population = 0

    def infectious_count(self):
        return self.infectious_nt + self.infectious_t + self.infectious_a

    def total_pop(self):
        return sum([self.susceptible, self.latent, self.infectious_nt,
                    self.infectious_t, self.infectious_a, self.recovered])


class SubpopulationNetwork(nx.DiGraph):
    """
    Network with subpopulation nodes and commuting edges between edges.
    Responsible for updating the effective properties of subpopulations.
    """

    def __init__(self):
        self.seasonal_variable = 1  # This is actually a constant, as we don't
                                    # account for seasonality.
        return

    def get_all_neighbors(self, node):
        """
        Returns all neighboring nodes, along with their associated data.
        """
        return [(n, self.node[n]['pop']) for n in nx.all_neighbors(self, node)]

    def set_commuting_rates(self, commuting_matrix):
        self.commuting_rates = commuting_matrix

    def set_commuting_return_rate(self, return_rate):
        self.return_rate = return_rate

    def add_subpopulation(self, subpop):
        self.add_node(subpop.id, pop=subpop)

    """
    To add an adge just use the nx.DiGraph method, like this:
    a = Subpopulation([params_for_a])
    b = Subpopulation([params_for_b])
    sn = SubpopulationNetwork()
    sn.add_subpopulation(a)
    sn.add_subpopulation(b)
    sn.add_edge(a.id, b.id)  <<<<
    """

    def get_exit_rate(self, node):
        """
        Return the commuting exit rate for a given node, as a function of
        commuting and return rates.
        """
        return (sum(self.commuting_rates[node]) / self.return_rate)

    def update_effective_populations(self):
        """
        For each subpopulation, update effective total population, taking into
        account the commuting rates between neighboring subpopulations.
        """
        for j, j_pop in self.nodes_iter(data=True):

            j_exit_rate = self.get_exit_rate(j)
            local_pop = (j_pop.infectious_nt +
                         (j_pop.total_pop - j_pop.infectious_nt) / (1 + j_exit_rate))

            other_pop = 0
            for i, i_pop in self.get_all_neighbors(j):
                i_exit_rate = sum(self.commuting_rates[i]) / self.return_rate
                other_pop += (
                              ((i_pop.total_pop - i_pop.infectious_nt) / 
                                (1 + i_exit_rate)
                              ) *
                              self.commuting_rates[i][j] / self.return_rate
                             )

            j_pop.effective_population = local_pop + other_pop

    def update_effective_forces_of_infection(self):
        """
        For each subpopulation, update the value of the effective force of
        infection parameter. 
        """
        for j, j_pop in self.nodes_iter(data=True):

            j_exit_rate = self.get_exit_rate(j)
            

class Model(object):
    """
    Contains model parameters and subpopulations network.
    Responsible for computing stochastic variables and running infection processes.
    """
    def __init__(self, subpop_network=None, p_exit_latent=None, p_recovery=None):
        if subpop_network is None:
            raise ValueError("Missing subpop_network argument")
        if p_exit_latent is None:
            raise ValueError("Missing p_exit_latent argument")
        if p_recovery is None:
            raise ValueError("Missing p_recovery argument")

        self.p_travel_allowed = self.generate_probability_travel_allowed()
        self.infection_transmission_rate = self.generate_infection_transmission_rate()
        self.p_asymptomatic_infection = self.generate_probability_asymptomatic()
        self.p_exit_latent = p_exit_latent
        self.p_recovery = p_recovery

        self.subpop_network = subpop_network

    # TODO: Implement this
    def generate_probability_asymptomatic(self):
        return 1

    # TODO: Implement this
    def generate_transmission_rate(self):
        return 



    def set_seasonality_variable(self, value):
        self.subpop_network.seasonal_variable = value


class CompartmentModel(object):

    def __init__(self):
        pass

    def compute_infection_rate(self):
        return (self.params['transmission_rate'] * self.infectious_count() /
                sum(self.comps.values()))

    def compute_commuting_rate

    def get_state(self):
        """
        Returns a copy of the model compartments.
        """
        return deepcopy(self.comps)

    def infectious_count(self):
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
        for x in ['susceptible', 'latent', 'symptomatic_no_travel',
                  'symptomatic_travel', 'asymptomatic', 'recovered']:
            if x not in comps.keys():
                raise Exception("GleamModel copartments must include {}".format(x))
        self.comps = comps

    def set_parameters(self, params):
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
