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
import pandas


def compute_infection_rate(comps, v):
    return v['transmission_rate'] * infectious_count(comps) / v['total_pop']


def infect(comps, v):
    """
    Updates counts in compartements (comps) according to model variables (v).
    """

    inf_rate = compute_infection_rate(comps, v)
    new_latent = comps['susceptible'] * inf_rate
    new_symptomatic_no_travel = (comps['latent'] *
                                 v['p_exit_latent'] *
                                 (1 - v['p_asymptomatic_infection']) *
                                 (1 - v['p_travel_permission']))
    new_symptomatic_travel = (comps['latent'] *
                              v['p_exit_latent'] *
                              (1 - v['p_asymptomatic_infection']) *
                              v['p_travel_permission'])
    new_asymptomatic = (comps['latent'] *
                        v['p_exit_latent'] *
                        v['p_asymptomatic_infection'])
    new_snt_recovered = comps['symptomatic_no_travel'] * v['p_recovery']
    new_st_recovered = comps['symptomatic_travel'] * v['p_recovery']
    new_asym_recovered = comps['asymptomatic'] * v['p_recovery']

    comps['susceptible'] -= new_latent
    comps['latent'] += new_latent
    comps['latent'] -= (new_symptomatic_no_travel + new_symptomatic_travel +
                        new_asymptomatic)
    comps['symptomatic_no_travel'] += new_symptomatic_no_travel
    comps['symptomatic_travel'] += new_symptomatic_travel
    comps['asymptomatic'] += new_asymptomatic
    comps['symptomatic_no_travel'] -= new_snt_recovered
    comps['symptomatic_travel'] -= new_st_recovered
    comps['asymptomatic'] -= new_asym_recovered
    comps['recovered'] += new_snt_recovered + new_st_recovered + new_asym_recovered


def infectious_count(comps):
    return (comps['symptomatic_travel'] +
            comps['symptomatic_no_travel'] +
            comps['asymptomatic'])


def print_graph(comps):
    """
    Prints a horizontal bar chart of the current compartment counts.
    """
    print("""
          susceptible:   {}  {}
          latent:        {}  {}
          sym_no_travel: {}  {}
          sym_travel:    {}  {}
          asymptomatic:  {}  {}
          recovered:     {}  {}
          """.format('-' * int((comps['susceptible'] // 13)),
                     int(comps['susceptible']),
                     '-' * int((comps['latent'] // 13)),
                     int(comps['latent']),
                     '-' * int((comps['symptomatic_no_travel'] // 13)),
                     int(comps['symptomatic_no_travel']),
                     '-' * int((comps['symptomatic_travel'] // 13)),
                     int(comps['symptomatic_travel']),
                     '-' * int((comps['asymptomatic'] // 13)),
                     int(comps['asymptomatic']),
                     '-' * int((comps['recovered'] // 13)),
                     int(comps['recovered']))
          )

model_variables = {'p_exit_latent': 0.8,
                   'transmission_rate': 0.5,
                   'total_pop': 2000,
                   'p_asymptomatic_infection': 0.1,
                   'p_travel_permission': 0.5,
                   'p_recovery': 0.2}

compartments = {'susceptible': 1950,
                'latent': 50,
                'symptomatic_no_travel': 0,
                'symptomatic_travel': 0,
                'asymptomatic': 0,
                'recovered': 0}


if __name__ == '__main__':
    print_graph(compartments)
    results = []
    for i in range(1000):
        results.append(deepcopy(compartments))
        # print('total pop: {}'.format(sum(compartments.values())))
        infect(compartments, model_variables)
    print_graph(compartments)
    with open('results.json', 'w') as f:
        json.dump(results, f)

