# !/usr/bin/env python
#  encoding: utf-8

import os
import pandas as pd
import pygal
import networkx as nx
import json
import re
import csv
import sys
from datetime import date, datetime
from collections import OrderedDict


def compute_commuting_flows(input_file, output_file):
    """
    Commuter flow estimation done using the radiation model of traffic flow.
    See article:
    Simini, Filippo, Marta C. González, Amos Maritan, and Albert-László Barabási. 2012.
    “A Universal Model for Mobility and Migration Patterns.”
    Nature 484 (7392): 96–100.
    """

    print('reading file')
    g = nx.read_graphml(input_file)

    # Percentage of commuters in total population
    commuter_percentage = 0.11
    counter = 1
    total_nodes = len(g.nodes())
    for i in g.nodes_iter():
        sys.stdout.write("\rProcessing node %d/%d" % (counter, total_nodes) )
        counter += 1
        sys.stdout.flush()
        pop_i = g.node[i]['pop']
        neighbors = set(nx.neighbors(g, i))
        for j in neighbors:
            pop_j = g.node[j]['pop']
            other_neighbors = neighbors - set([j])
            radius = g.edge[i][j]['Total_Length']
            others_in_radius = [nb for nb in other_neighbors
                                if g.edge[i][nb]['Total_Length'] < radius]
            pop_in_radius = sum([g.node[o]['pop'] for o in others_in_radius])
            g.edge[i][j]['commuting_rate'] = (
                                              (pop_i * pop_j) /
                                              (
                                                (pop_i + pop_in_radius) * 
                                                (pop_i + pop_j + pop_in_radius)
                                              )
                                             ) * commuter_percentage
    print('\nwriting to file')
    nx.write_graphml(g, output_file)


def compute_effective_population(input_file, output_file):
    """
    For subpopulation node of a graph, compute effective total population,
    taking into account the commuting rates between neighboring subpopulations.
    """

    g = nx.read_graphml(input_file)

    for node_id in g.nodes_iter(data=True):
        local_pop = g.node[node_id]['pop'] / (1 + g.node[node_id]['sigma_by_tau'])
        other_pop = sum([g.node[neighbor]['pop'] *
                           (g.edge[neighbor][node_id]['sigma_prop_by_tau'] /
                            (1 + g.node[neighbor]['sigma_by_tau'])
                           )
                         for neighbor in g.predecessors(node_id)])
        g.node[node_id]['effective_population'] = local_pop + other_pop

    nx.write_graphml(g, output_file)


def compute_sigmas(input_file, output_file, tau):
    """
    Add sigma_by_tau and sigma_prop_by_tau attributes to vertices of a graph.
    """
    g = nx.read_graphml(input_file)

    for source in g.nodes_iter():
        sigma_by_tau = sum([g.edge[source][dest]['commuting_rate']
                            for dest in g.successors(source)]) / tau
        for dest in g.successors(source):
            g.edge[dest]['sigma_prop_by_tau'] = g.edge[source][dest]['commuting_rate'] / tau

    nx.write_graphml(g, output_file)


def format_graph(input_graph_filepath, output_graph_filepath):
    """
    This is what I used to format the rwa_net_pruned.graphml into what pygleam
    expects.
    """
    g = nx.read_graphml(input_graph_filepath)
    for n in g.nodes():
        g.node[n] = {'compartments': {'susceptible': 0,
                                      'latent': 0,
                                      'infectious_a': 0,
                                      'infectious_t': 0,
                                      'infectious_nt': 0,
                                      'recovered': 0},
                     'coordinates': (g.node[n]['lat'], g.node[n]['lon']),
                     'name': g.node[n]['name']}
    nx.write_gpickle(g, output_graph_filepath)


def generate_geojson_base_nodes(input_file, output_file):
    """
    Given a network file, will collect latitude and longitude for every node and
    generate a geojson file to allow mapping of the nodes on a Leaflet map.
    """
    def generate_point_object(node):
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
                       "population": node['pop']
                   }
               }

    g = nx.read_graphml(input_file)
    output = [generate_point_object(node[1]) for node in g.nodes_iter(data=True)]
    output_str = 'base_nodes = ' + json.dumps(output)

    with open(output_file, 'w') as f:
        f.write(output_str)


def generate_geojson_commuting_edges(input_file, output_file, min_commuting_rate):
    """
    Given a network file, will extract all edges with commuting rate >= to
    commuting_rate, and will generate a geojson file to allow mapping of the
    edges on a Leaflet map.
    """
    def generate_line_object(coords_1, coords_2, comm_rate):
        return {
                "type": "Feature",
                "geometry": {
                              "type": "LineString",
                              "coordinates": [
                                                [coords_1[0], coords_1[1]],
                                                [coords_2[0], coords_2[1]]
                                             ]
                            },
                "properties": {
                                "commuting_rate": comm_rate
                              }
               }
    g = nx.read_graphml(input_file)
    edges = [x for x in g.edges(data=True) if x[2]['commuting_rate'] >= 0.05]
    output = []
    for e in edges:
        coords_1 = [g.node[e[0]]['lon'], g.node[e[0]]['lat']]
        coords_2 = [g.node[e[1]]['lon'], g.node[e[1]]['lat']]
        comm_rate = e[2]['commuting_rate']
        output.append(generate_line_object(coords_1, coords_2, comm_rate))
    output_str = 'edges = ' + json.dumps(output)

    with open(output_file, 'w') as f:
        f.write(output_str)


def get_csv_data_from_results_by_node(results_filepath, output_csv_filepath):
    """
    Takes a GeoJSON result file and outputs a CSV file with following row structure:
    node_name, timestep, susceptibles, latents, infectious, recovered
    """

    with open(results_filepath) as f:
        raw_contents = f.read()
        json_str = re.sub(r'\w+ = ', '',raw_contents)
        data = json.loads(json_str)
        nodes_dict = {}
        for obj in data:
            name = int(obj['properties']['name'])
            times = obj['properties']['times']
            compartments = obj['properties']['compartments']
            if name not in nodes_dict.keys():
                nodes_dict[name] = {}
            for t in times:
                date_obj = datetime.strptime(t, '%Y-%m-%d')
                nodes_dict[name][date_obj] = compartments
        
        one_key = ''
        for key in nodes_dict.keys():
            nodes_dict[key] = OrderedDict(sorted(nodes_dict[key].items()))
            one_key = key

        nodes_dict = OrderedDict(sorted(nodes_dict.items()))

        timesteps = range(1, len(nodes_dict[one_key])+1)
        
    with open(output_csv_filepath, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['node_name', 'timestep', 'susceptible',
                         'latent', 'infectious', 'recovered'])
        for name in nodes_dict.keys():
            timestep = 1
            for date in nodes_dict[name].keys():
                row = [name, timestep,
                       nodes_dict[name][date]['susceptible'],
                       nodes_dict[name][date]['latent'],
                       (nodes_dict[name][date]['infectious_a'] +
                        nodes_dict[name][date]['infectious_t'] +
                        nodes_dict[name][date]['infectious_nt']),
                       nodes_dict[name][date]['recovered']]
                writer.writerow(row)
                timestep += 1


def get_csv_data_from_results_global(results_filepath, output_csv_filepath):
    """
    Takes a GeoJSON result file and outputs a CSV file with following row structure:
    timestep, susceptibles, latents, infectious, recovered
    """
    data = get_global_compartment_values_by_timestep(results_filepath)
    with open(output_csv_filepath, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['timestep', 'susceptible', 'latent', 'infectious', 'recovered'])
        for i in range(len(data['timesteps'])):
            writer.writerow([data['timesteps'][i],
                             data['sus'][i],
                             data['lat'][i],
                             data['inf'][i],
                             data['rec'][i]]
                            )


# TODO: convert result file to CSV file with following row structure:
# node_name, timestep, susceptibles, latents, infectious_a, infectious_nt, infectious_t, recovered

def get_global_compartment_values_by_timestep(results_filepath):
    """
    From a GeoJSON results file, get total susceptible, infected, latent and recovered
    for each timestep.
    """

    with open(results_filepath) as f:
        raw_contents = f.read()
        json_str = re.sub(r'\w+ = ', '',raw_contents)
        data = json.loads(json_str)
        date_comps = {}
        for obj in data:
            times = obj['properties']['times']
            compartments = obj['properties']['compartments']
            for t in times:
                date_obj = datetime.strptime(t, '%Y-%m-%d')
                if date_obj not in date_comps.keys():
                    date_comps[date_obj] = {
                                             "infectious_t": 0,
                                             "infectious_a": 0,
                                             "recovered": 0,
                                             "infectious_nt": 0,
                                             "latent": 0,
                                             "susceptible": 0
                                           }
                for c in date_comps[date_obj].keys():
                    date_comps[date_obj][c] += compartments[c]

        s = sorted(date_comps.items())
        timesteps = range(1, len(s)+1)
        sus = [x[1]['susceptible'] for x in s]
        lat = [x[1]['latent'] for x in s]
        inf = [x[1]['infectious_a'] + x[1]['infectious_t'] + x[1]['infectious_nt']
               for x in s]
        rec = [x[1]['recovered'] for x in s]

        outbreak_end = rec.index(max(rec))
        for x in [timesteps, sus, lat, inf, rec]:
            x = x[:outbreak_end+1]

        return {'sus': sus, 'lat': lat, 'inf': inf, 'rec': rec, 'timesteps': timesteps}


def get_recovered_counts_from_results(input_folder):
    """
    Retrieve recovered counts at last time steps from many simulations,
    possibly contained in many files in the input_folder.
    """
    counts = []
    sus_counts = []
    for file in os.listdir(input_folder):
        data = json.load(open(os.path.join(input_folder, file)))
        for simul in data:
            total_rec = 0
            total_sus = 0
            for history in simul.values():
                total_sus += history['susceptible'][0]
                # Take the number of recovered at last timestep
                total_rec += history['recovered'][-1]
            counts.append(total_rec)
            sus_counts.append(total_sus)
    print(counts)
    print(sus_counts)


def graphml_to_geojson(graphml_input_filepath, geojson_output_filepath):
    """
    Takes a graphml file as input and writes a correctly formatted geojson
    file.
    """

    def geojson_builder(node):
        """
        Takes a node from the graph as input and returns the data in valid
        GeojSON format.
        """
        return {
                'type': 'Feature',
                'properties': {
                                'name': node['name'],
                                'population': node['pop']
                              },
                'geometry': {
                                'type': 'Point',
                                'coordinates': [node['lat'], node['lon']]
                            }
                }

    G = nx.read_graphml(graphml_input_filepath)

    geojson_list = [geojson_builder(G.node[x]) for x in G.nodes_iter()]

    file_content = 'nodes = ' + json.dumps(geojson_list)

    with open(geojson_output_filepath, 'w') as f:
        f.write(file_content)


def plot_epidemic_curve_from_results(results_filepath, plot_output_filepath):
    """
    Takes a GeoJSON result file and plot the global epidemic curve from it.
    """
    
    comps = get_global_compartment_values_by_timestep(results_filepath)

    custom_style = pygal.style.DarkStyle
    custom_style.colors = ["#feed6c", "#bf4646", "#56c2d6", "#516083"]
    line_chart = pygal.Line(dots_size=1,
                            show_dots=False,
                            show_y_guides=False,
                            legend_at_bottom_columns=2,
                            x_labels_major_every=5,
                            x_label_rotation=270,
                            show_minor_x_labels=False,
                            truncate_label=-1,
                            margin_bottom=50,
                            style=custom_style)

    line_chart.title = 'Simulation results'
    line_chart.x_labels = dates
    line_chart.add('susceptible', comps['sus'])
    line_chart.add('latent', comps['lat'])
    line_chart.add('infectious', comps['inf'])
    line_chart.add('recovered', comps['rec'])
                                   
    line_chart.render_in_browser()
    line_chart.render_to_file(plot_output_filepath)


def plot_histogram(input_file):
    style.use('dark_background')
    with open(input_file) as f:
        counts = json.load(f)

    minimum, maximum = min(counts), max(counts)

    fig, ax = pl.subplots(1,1, tight_layout=True)
    counts, bins, patches = ax.hist(counts, bins=10 ** np.linspace(np.log10(minimum),
                                                                   np.log10(maximum), 70),
                                    facecolor='green',
                                    edgecolor='black')

    patches[-1].set_color('red')

    # pl.figure()
    # h = pl.hist(counts, bins = 10 ** np.linspace(np.log10(minimum), np.log10(maximum), 70))
    pl.gca().set_xscale("log")
    ax.set_xlabel('Total number of infected during outbreak', fontsize=15)
    ax.set_ylabel('Number of simulations', fontsize=15)
    ax.tick_params(labelsize=10)
    ax.set_ylim(0, 40)

    pl.savefig(input_file[-14:-6]+'.png')


def prune_edges_with_min_cr(input_filepath, output_filepath, min_cr):
    """
    Removes from graph the edges that have a commuting rate under min_cr.
    """
    print('loading graph')
    graph = nx.read_graphml(input_filepath)
    print('Identifying bad edges')
    bad_edges = [x for x in graph.edges_iter(data=True)
                 if x[2]['commuting_rate'] < min_cr]

    print('removing bad edges')
    graph.remove_edges_from(bad_edges)
    print('writing results to disk')
    nx.write_graphml(graph, output_filepath)


if __name__ == '__main__':
    # G = nx.read_graphml('/data/influenza/rwanda/rwa_net.graphml')
    # prune_edges_with_max_distance(G, 50000)
    # nx.write_graphml(G, '/data/influenza/rwanda/rwa_net_pruned.graphml')
    
    input_file = '/home/hugo/data/pygleam/2016-07-15_rwa-net.graphml'
    output_file = '/home/hugo/data/pygleam/rwa-net_cr_wrong.graphml'
    compute_commuting_flows(input_file, output_file)

    input_file = '/home/hugo/data/pygleam/rwa-net_cr_wrong.graphml'
    output_file = '/home/hugo/data/pygleam/rwa-net_cr_pruned_wrong.graphml'
    prune_edges_with_min_cr(input_file, output_file, 0.001)

    # input_file = '/home/hugo/Projects/gleam/data/rwa_net.graphml'
    # output_file = '/home/hugo/Projects/gleam/data/base_nodes.jsonp'
    # generate_geojson_base_nodes(input_file, output_file)

    # plot_epidemic_curve_from_results('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
    #                                  'charts/H1N1_node-890_seed-1_n-100_vacc-0.0.svg')

    # get_csv_data_from_results_global('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
    #                                  'csv/H1N1_node-890_seed-1_n-100_vacc-0.0_GLOBAL.csv')

    # get_csv_data_from_results_by_node('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
    #                                   'csv/H1N1_node-890_seed-1_n-100_vacc-0.0_BYNODE.csv')

    # get_recovered_counts_from_results('output/test')

    # plot_histogram('/home/hugo/Projects/pygleam/output/recovered_counts_0.0.json')
    # plot_histogram('/home/hugo/Projects/pygleam/output/recovered_counts_0.2.json')
    # plot_histogram('/home/hugo/Projects/pygleam/output/recovered_counts_0.4.json')
    # plot_histogram('/home/hugo/Projects/pygleam/output/recovered_counts_0.6.json')
    # plot_histogram('/home/hugo/Projects/pygleam/output/recovered_counts_0.8.json')