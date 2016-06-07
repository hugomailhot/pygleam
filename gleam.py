# !/usr/bin/env python
# encoding: utf-8
"""
Implementation of the GLEaM model. See #reference article:

Balcan, Duygu, Bruno Gonçalves, Hao Hu, José J. Ramasco, Vittoria Colizza,
and Alessandro Vespignani. 2010.
“Modeling the Spatial Spread of Infectious Diseases: The GLobal Epidemic and
Mobility Computational Model.”
Journal of Computational Science 1 (3): 132–45.
"""

from collections import Counter
from copy import deepcopy
import math
import numpy as np
import networkx as nx
from datetime import date, timedelta
import json
from time import strftime
from pprint import pprint

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
        Args:
            subpop_network (nx.DiGraph): graph representation of the subpopulations
                                         and their commuting relationships
            params (dict): model parameters

        Raises:
            ValueError: Description
            ValueError
            param dict does not contain all required values
        """
        for param in ['p_exit_latent', 'p_recovery', 'p_asymptomatic',
                      'p_travel_allowed', 'commuting_return_rate', 'asym_downscaler']:
            if param not in params.keys():
                raise ValueError("Missing {} parameter".format(param))

        super(Model, self).__init__(subpop_network)

        self.p_exit_latent = params['p_exit_latent']
        self.p_recovery = params['p_recovery']
        self.p_asymptomatic = params['p_asymptomatic']
        self.p_travel_allowed = params['p_travel_allowed']
        self.return_rate = params['commuting_return_rate']
        self.asym_downscaler = params['asym_downscaler']
        self.starting_date = params['starting_date']


        for i in self.nodes_iter():
            self.node[i]['pop'] = math.ceil(self.node[i]['pop'])
            self.node[i]['compartments'] = Counter(
                                           {'susceptible': math.ceil(self.node[i]['pop']),
                                            'latent': 0,
                                            'infectious_a': 0,
                                            'infectious_t': 0,
                                            'infectious_nt': 0,
                                            'recovered': 0
                                            })
            # Store state at beginning.
            self.node[i]['history'] = [deepcopy(self.node[i]['compartments'])]
            self.node[i]['exit_rate'] = self.get_exit_rate(i)

        for i in self.nodes_iter():
            self.node[i]['effective_population'] = self.effective_population(i)

    def compute_long_distance_travels(self):
        """
        Redistributes travelers among neighboring nodes, according to
        long-distance traffic information.
        """
        # TODO: implement this
        pass

    def draw_new_latent_count(self, node_id):
        """
        For a given node, draws a random value from a binomial distribution
        defined by susceptible population and chance of contagion.

        Args:
            node_id (int): id of the nx.DiGraph node

        Returns:
            float: Probability extracted from the binomial distribution
        """
        if self.node[node_id]['compartments']['susceptible'] < 1:
            return 0
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
        if self.node[node_id]['compartments']['latent'] < 1:
            return (0, 0, 0)
        p_a = self.p_asymptomatic
        p_t = (1 - self.p_asymptomatic) * self.p_travel_allowed
        p_nt = (1 - self.p_asymptomatic) * (1 - self.p_travel_allowed)
        a, t, nt = np.random.multinomial(self.node[node_id]['compartments']['latent'],
                                         [p_a, p_t, p_nt], size=1)[0]
        return [math.ceil(x * self.p_exit_latent) for x in [a, t, nt]]

    def draw_new_recovered_counts(self, node_id):
        """
        For a given node, draw from a binomial distribution the number of infectious
        people that will recover in next time step.
        """
        node_pop = self.node[node_id]['compartments']
        if (node_pop['infectious_a'] +
            node_pop['infectious_t'] +
                node_pop['infectious_nt']) < 1:
            return (0, 0, 0)
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
        node_exit_rate = self.node[node_id]['exit_rate']
        local_pop = (node_pop['infectious_nt'] +
                     (sum(node_pop.values()) - node_pop['infectious_nt']) /
                     (1 + node_exit_rate))

        other_pop = 0
        # Only consider neighbors that have an edge inbound to local node, ie predecessors
        for nb_id in self.predecessors(node_id):
            nb_pop = self.node[nb_id]['compartments']
            nb_exit_rate = self.node[nb_id]['exit_rate']
            other_pop += (((sum(nb_pop.values()) - node_pop['infectious_nt']) /
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
        local_exit_rate = self.node[node_id]['exit_rate']
        local_foi = self.force_of_infection(node_id) / (1 + local_exit_rate)

        nbs_foi = 0
        # Only consider neighbors that can be attained from local node, ie successors
        for nb_id in self.successors(node_id):
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
        local_exit_rate = self.node[node_id]['exit_rate']
        local_infectious = (node_pop['infectious_nt'] +
                            (node_pop['infectious_t'] + scaled_asym_pop) /
                            (1 + local_exit_rate))

        neighbors_infectious = 0
        # Only consider neighbors that have an edge inbound to local node, ie predecessors
        for nb_id in self.predecessors(node_id):
            nb_pop = self.node[nb_id]['compartments']
            scaled_asym_pop = self.asym_downscaler * nb_pop['infectious_a']
            nb_exit_rate = self.node[nb_id]['exit_rate']
            commuting_rate = self.edge[nb_id][node_id]['commuting_rate']
            commuting_nb_inf = ((nb_pop['infectious_t'] + scaled_asym_pop) /
                                (1 + nb_exit_rate)) * commuting_rate
            neighbors_infectious += commuting_nb_inf

        total_infectious = local_infectious + neighbors_infectious
        return (self.seasonality() / self.node[node_id]['effective_population'] *
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

            compartments = self.node[node_id]['compartments']

            compartments['susceptible'] -= new_latent
            compartments['latent'] += new_latent
            compartments['latent'] -= total_new_inf
            compartments['infectious_a'] += new_inf_a
            compartments['infectious_a'] -= new_a_to_r
            compartments['infectious_t'] += new_inf_t
            compartments['infectious_t'] -= new_t_to_r
            compartments['infectious_nt'] += new_inf_nt
            compartments['infectious_nt'] -= new_nt_to_r
            compartments['recovered'] += total_new_recovered

            self.node[node_id]['history'].append(Counter(compartments))

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
        return 1.5 * self.p_recovery

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

    def run_n_simulations(self, n, timesteps=200):
        """Using the same starting conditons, will run the infection process
        n times, then for each node will place in history the average value over
        all simulations for each time step.

        Args:
            n (int): Number of simulations to run
        """
        def there_is_infected_nodes(model):
            return any(any(model.node[node_id]['compartments'][comp] != 0
                           for comp in ['latent', 'infectious_t',
                                       'infectious_a', 'infectious_nt'])
                       for node_id in model.nodes_iter())
        self.fresh_copy = deepcopy(self)
        new_compartment = {'susceptible': 0,
                           'latent': 0,
                           'infectious_a': 0,
                           'infectious_t': 0,
                           'infectious_nt': 0,
                           'recovered': 0}

        for node_id in self.nodes_iter():
            self.node[node_id]['history'] = []
            for i in range(timesteps):
                self.node[node_id]['history'].append(Counter(new_compartment))

        for i in range(1, n + 1):
            new_model = deepcopy(self.fresh_copy)
            print(strftime('%H:%M:%S') + '  Simulation #{}'.format(i))
            for t in range(timesteps):
                new_model.infect()
            # Add compartment values for each node, divided by number of iterations,
            # to model node histories
            for node_id in new_model.nodes_iter():
                history = new_model.node[node_id]['history']
                for i in range(timesteps):
                    weighted_comps = Counter({k: v / n for k, v in history[i].items()})
                    self.node[node_id]['history'][i] += weighted_comps

    def generate_timestamped_geojson_output(self, output_file):
        """Generates the file to be read by the Leaflet.js visualization script.

        Args:
            output_file (str): Path of the output file
        """

        def generate_point_object(node, state, date):
            return {
                        "type": "Feature",
                        "geometry": {
                            "type": "Point",
                            "coordinates": [
                                node['lon'],
                                node['lat']
                            ]
                        },
                        "properties": {
                            "name": node['name'],
                            "population": node['pop'],
                            "times": [str(date)],
                            "compartments": {
                                "susceptible": round(state['susceptible']),
                                "latent": round(state['latent']),
                                "infectious_t": round(state['infectious_t']),
                                "infectious_nt": round(state['infectious_nt']),
                                "infectious_a": round(state['infectious_a']),
                                "recovered": round(state['recovered'])
                            }
                        }
                    }

        output = []
        for n_id in self.nodes_iter():
            node = self.node[n_id]
            history = self.node[n_id]['history']
            date = self.starting_date
            output.append(generate_point_object(node, history[0], date))

            for i in range(1, len(history) - 1):
                date += timedelta(days=1)
                if history[i] == history[i-1]:
                    output[-1]['properties']['times'].append(str(date))
                else:
                    output.append(generate_point_object(node, history[i], date))

        output_str = 'nodes = ' + json.dumps(output)
        with open(output_file, 'w') as f:
            f.write(output_str)


if __name__ == '__main__':

    # Define model parameters
    starting_date = date(2016, 7, 11)
    model_parameters = {'p_exit_latent': 1 / 1.1,
                        'p_recovery': 1 / 2.95,
                        'p_asymptomatic': 0.2,
                        'p_travel_allowed': 0.3,
                        'commuting_return_rate': 3,
                        'asym_downscaler': 1,
                        'starting_date': starting_date}
    
    graph_filepath = 'data/rwa_net.graphml'
    gleam = Model(nx.read_graphml(graph_filepath), model_parameters)

    # Kigali is n842
    starting_node = 'n190'
    seeds = 5
    gleam.seed_infectious(starting_node, seeds=seeds)
    
    nb_simulations = 20
    timesteps_per_simul = 200
    gleam.run_n_simulations(n=nb_simulations, timesteps=timesteps_per_simul)

    output_file = 'output/node-{}_seed-{}_n-{}.jsonp'.format(starting_node[1:],
                                                             seeds,
                                                             nb_simulations)
    gleam.generate_timestamped_geojson_output(output_file)
