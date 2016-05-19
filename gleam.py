"""
Implementationg of the GLEaM model. See reference article:

Balcan, Duygu, Bruno Gonçalves, Hao Hu, José J. Ramasco, Vittoria Colizza,
and Alessandro Vespignani. 2010.
“Modeling the Spatial Spread of Infectious Diseases: The GLobal Epidemic and
Mobility Computational Model.”
Journal of Computational Science 1 (3): 132–45.
"""
# !/usr/bin/env python
# encoding: utf-8


from copy import deepcopy
import math
import pandas as pd
import pygal
from pprint import pprint
import random
from utilities import plot_results
import numpy as np
import networkx as nx


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
    # TODO: revise every exit_rate usage. Ensure 1 is added to it each time.

    def __init__(self, subpop_network, params):
        """
        Args:
            subpop_network (nx.DiGraph): graph representation of the subpopulations 
                                         and their commuting relationships
            params (dict): model parameters
        
        Raises:
            ValueError: Description
            ValueError
            param dict does not contain all required values
        """
        super(Model, self).__init__(subpop_network)

        for param in ['p_exit_latent', 'p_recovery', 'p_asymptomatic',
                      'p_travel_allowed', 'commuting_return_rate', 'asym_downscaler']:
            if param not in params.keys():
                raise ValueError("Missing {} parameter".format(param))

        self.p_exit_latent = params['p_exit_latent']
        self.p_recovery = params['p_recovery']
        self.p_asymptomatic = params['p_asymptomatic']
        self.p_travel_allowed = params['p_travel_allowed']
        self.commuting_return_rate = params['commuting_return_rate']
        self.asym_downscaler = params['asym_downscaler']

    def draw_new_latent_count(self, node_id):
        """
        For a given node, draws a random value from a binomial distribution
        defined by susceptible population and chance of contagion.
        
        Args:
            node_id (int): id of the nx.DiGraph node
        
        Returns:
            float: Probability extracted from the binomial distribution
        """
        return np.random.binomial(self.node[node_id]['compartments']['susceptible'],
                                  self.effective_force_of_infection(node_id))

    def draw_new_infectious_counts(self, node_id):
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
        a, t, nt = np.random.multinomial(self.node[node_id]['compartments']['latent'],
                                         [p_a, p_t, p_nt], size=1)
        return [x * self.p_exit_latent for x in [a, t, nt]]

    def draw_new_recovered_counts(self, node_id):
        """
        For a given node, draw from a binomial distribution the number of infectious
        people that will recover in next time step.
        """
        node_pop = self.node[node_id]['compartments']
        a_recovered = np.random.binomial(node_pop['infectious_a'], self.p_recovery)
        t_recovered = np.random.binomial(node_pop['infectious_t'], self.p_recovery)
        nt_recovered = np.random.binomial(node_pop['infectious_nt'], self.p_recovery)
        return (a_recovered, t_recovered, nt_recovered)

    def effective_population(self, node_id):
        """
        For given subpopulation node, compute effective total population,
        taking into account the commuting rates between neighboring subpopulations.
        
        Args:
            node_id (int): id of the nx.DiGraph node
        
        Returns:
            int: Effective population of given node.
        """
        node_pop = self.node[node_id]['compartments']
        node_exit_rate = self.get_exit_rate(node_id)
        local_pop = (node_pop['infectious_nt'] +
                     (sum(node_pop.values()) - node_pop['infectious_nt']) /
                     (1 + node_exit_rate))

        other_pop = 0
        for nb_id in nx.neighbors(self, node_id):
            nb_pop = self.node[nb_id]['compartments']
            nb_exit_rate = self.get_exit_rate(nb_id)
            other_pop += (((sum(nb_pop.values()) - nb_pop.infectious_nt) /
                           (1 + nb_exit_rate)) *
                          self.edge[nb_id][node_id]['commuting_rate'] / self.return_rate)

        return local_pop + other_pop

    def effective_force_of_infection(self, node_id):
        """For given subpopulation node, computes effective for of infection,
        
        taking into account two terms:
            1)  the local force of infection;
            2)  the forces of infection in neighboring nodes, scaled respectively
                by the amount of people from local node that commute to these neighbors.
        
        Args:
            node_id (int): id of the nx.DiGraph node
        
        Returns:
            float: Effective force of infection at given node.
        """
        local_exit_rate = self.get_exit_rate(node_id)
        local_foi = self.force_of_infection(node_id) / (1 + local_exit_rate)

        nbs_foi = 0
        for nb_id in nx.neighbors(self, node_id):
            commuting_rate_local_to_nb = self.edge[node_id][nb_id]['commuting_rate']
            effective_commuting_rate = commuting_rate_local_to_nb / self.return_rate
            nbs_foi += (self.force_of_infection(nb_id) * effective_commuting_rate /
                        (1 + local_exit_rate))

        return local_foi + nbs_foi

    def force_of_infection(self, node_id):
        """For a given node, computes the force of infection, taking into two
        
        terms:
            1)  the number of infectious people in local node;
            2)  the number of infectious people from neighboring node that commute
                to local node.
        
        Args:
            node_id (int): id of the nx.DiGraph node
        
        Returns:
            float: Force of infection at given node.
        """
        node_pop = self.node[node_id]['compartments']
        scaled_asym_pop = self.asym_downscaler * node_pop['infectious_a']
        local_exit_rate = self.get_exit_rate(node_id)
        local_infectious = (node_pop['infectious_nt'] +
                            (node_pop['infectious_t'] + scaled_asym_pop) /
                            (1 + local_exit_rate))

        neighbors_infectious = 0
        for nb_id in nx.neighbors(self, node_id):
            nb_pop = self.node[nb_id]['compartments']
            scaled_asym_pop = self.asym_downscaler * nb_pop['infectious_a']
            nb_exit_rate = self.get_exit_rate(nb_id)
            commuting_rate = self.edge[nb_id][node_id]['commuting_rate']
            commuting_nb_inf = ((nb_pop['infectious_t'] + scaled_asym_pop) /
                                (1 + nb_exit_rate)) * commuting_rate
            neighbors_infectious += commuting_nb_inf

        total_infectious = local_infectious + neighbors_infectious
        return (self.seasonality() / self.effective_population(node_id) *
                total_infectious)

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

    def infect(self):
        """
        Runs one step of the infection process on every node in the network.
        """
        for node_id in self.nodes_iter():
            new_latent = self.draw_new_latent_count(node_id)
            new_inf_a, new_inf_t, new_inf_nt = self.draw_new_infectious_counts(node_id)
            new_a_to_r, new_t_to_r, new_nt_to_r = self.draw_new_recovered_counts(node_id)
            total_new_inf = new_inf_a + new_inf_t + new_inf_nt
            total_new_recovered = new_a_to_r + new_t_to_r + new_nt_to_r

            self.node[node_id]['compartments']['susceptible'] -= new_latent
            self.node[node_id]['compartments']['latent'] += new_latent
            self.node[node_id]['compartments']['latent'] -= total_new_inf
            self.node[node_id]['compartments']['infectious_a'] += new_inf_a
            self.node[node_id]['compartments']['infectious_a'] -= new_a_to_r
            self.node[node_id]['compartments']['infectious_t'] += new_inf_t
            self.node[node_id]['compartments']['infectious_t'] -= new_t_to_r
            self.node[node_id]['compartments']['infectious_nt'] += new_inf_nt
            self.node[node_id]['compartments']['infectious_nt'] -= new_nt_to_r
            self.node[node_id]['compartments']['recovered'] += total_new_recovered

    def compute_long_distance_travels(self):
        """
        Redistributes travelers among neighboring nodes, according to
        long-distance traffic information.
        """
        # TODO: implement this
        pass

    def seasonality(self, hemisphere=None):
        """
        Computes the scalar factor to apply on force of infection.
        
        Args:
            hemisphere (str): either 'north', 'south' or None
        
        Returns:
            float: seasonality value for the given hemisphere, or 1 if no hemisphere
                   value was given
        
        Raises:
            ValueError: Description
        """
        # TODO: implement seasonality computation
        if hemisphere not in ['north', 'south', None]:
            raise ValueError("hemisphere argument must be in ['north', 'south', None]")
        # if hemisphere is None:
        return 1

    def seed_infectious(self, node_id, seeds=1, inf_type='infectious_t'):
        """
        Transfers a number of people from susceptible compartment to an infectious
        compartment at given node_id.
        
        Args:
            node_id (str): Node in which we want to seed infectious people
            seeds (int): Number of people to be seeded. 1 by default.
            inf_type (str): Specifies the type of infectious to be seeded.
        """
        self.node[node_id]['compartments']['susceptible'] -= seeds
        self.node[node_id]['compartments'][inf_type] += seeds

    def total_infectious(self, node_id):
        """Returns the total number of infectious people in a given node population.
        
        Args:
            node_pop (dict): Contains population values for the different compartments
        
        Returns:
            int: Sum of values in infectious compartments
        """
        node_pop = self.node[node_id]['compartments']
        return sum([node_pop['infectious_nt'], node_pop['infectious_t'],
                    node_pop['infectious_a']])

    def total_pop(self, node_pop):
        """Returns the total population in a given node population.
        
        Args:
            node_pop (dict): Contains population values for the different compartments
        
        Returns:
            int: Sum of values in all compartments
        """
        return sum(node_pop.values())


if __name__ == '__main__':
    # parameters = {'p_exit_latent': 0.1,
    #               'transmission_rate': 0.5,
    #               'p_asymptomatic_infection': 0.1,
    #               'p_travel_permission': 0.3,
    #               'p_recovery': 0.05}

    # populations = {'susceptible': 1950,
    #                'latent': 0,
    #                'symptomatic_no_travel': 10,
    #                'symptomatic_travel': 0,
    #                'asymptomatic': 0,
    #                'recovered': 0}
    # gm = CompartmentModel()
    # gm.set_populations(populations)
    # gm.set_parameters(parameters)
    # results = []
    # for i in range(1000):
    #     results.append(gm.get_state())
    #     gm.infect()
    # save_file = '/home/hugo/Projects/networks/network_theory_application_course_project/visualization/svg_test.svg'
    # print
    # plot_results(results)
    model_parameters = {'p_exit_latent': 0.2,
                        'p_recovery': 1 / 6,
                        'p_asymptomatic': 0.2,
                        'p_travel_allowed': 0.3,
                        'commuting_return_rate': 3,
                        'asym_downscaler': 1}
    graph_filepath = '/data/influenza/rwanda/rwa_net_formatted.gpickle'
    gleam = Model(nx.read_gpickle(graph_filepath), model_parameters)
    gleam.seed_infectious('n587')
    for i in range(100):
        gleam.infect()
    print(gleam.node['n587']['compartments'])
    # TODO: update network with commuting rates