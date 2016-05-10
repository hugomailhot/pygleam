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
import json
import math


class GleamModel(object):

    def __init__(self):
        pass

    def compute_infection_rate(self):
        return (self.params['transmission_rate'] * self.infectious_count() /
                sum(self.comps.values()))

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

        inf_rate = self.compute_infection_rate(self.comps, self.params)
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

    def set_compartments(self, comps):
        for x in ['susceptible', 'latent', 'symptomatic_no_travel',
                  'symptomatic_travel', 'asymptomatic', 'recovered']:
            if x not in comps.keys():
                raise Exception("GleamModel copartments must include {}".format(x))
        self.comps = comps

    def set_model_params(self, params):
        for x in ['p_exit_latent', 'transmission_rate', 'p_asymptomatic_infection',
                  'p_travel_permission', 'p_recovery']:
            if x not in params.keys():
                raise Exception("GleamModel parameters must include {}".format(x))
        self.params = params







if __name__ == '__main__':
    parameters = {'p_exit_latent': 0.2,
                  'transmission_rate': 0.6,
                  'total_pop': 2000,
                  'p_asymptomatic_infection': 0.1,
                  'p_travel_permission': 0.5,
                  'p_recovery': 0.1}

    compartments = {'susceptible': 1950,
                    'latent': 0,
                    'symptomatic_no_travel': 10,
                    'symptomatic_travel': 0,
                    'asymptomatic': 0,
                    'recovered': 0}
    gm = GleamModel()
    gm.set_compartments(compartments)
    gm.set_parameters(parameters)
    results = []
    for i in range(1000):
        results.append(gm.get_state())
        # print('total pop: {}'.format(sum(compartments.values())))
        gm.infect()
    with open('results.json', 'w') as f:
        json.dump(results, f)

