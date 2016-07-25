# !/usr/bin/env python
#  encoding: utf-8

import pandas as pd
import pygal
import networkx as nx
import json
import re
import csv
from datetime import date, datetime
from collections import OrderedDict
# import folium
# from folium import plugins

# def plot_results(res, svg_filepath):
#     results = pd.DataFrame(res)
#     results['infectious'] = results[['infectious_nt',
#                                      'infectious_t',
#                                      'infectious_a']].sum(axis=1)

#     outbreak_peak = results['infectious'].idxmax()
#     outbreak_end = results['infectious'][outbreak_peak:].idxmin()
#     peak_inf_percentage = results['infectious'].max() / 2000.0

#     textstr = ('outbreak peak: {}\n'.format(outbreak_peak) +
#                'highest infection %%: {}%%\n'.format(peak_inf_percentage) +
#                'outbreak end: {}'.format(outbreak_end))

#     results = results[:outbreak_end]

#     custom_style = pygal.style.Style(colors=["#fdae61", "#d7191c", "#a6d96a", "#1a9641"])
#     custom_style = pygal.style.DarkStyle
#     custom_style.colors = ["#feed6c", "#bf4646", "#56c2d6", "#516083"]
#     line_chart = pygal.StackedLine(dots_size=1,
#                                    fill=True,
#                                    show_dots=False,
#                                    show_y_guides=False,
#                                    legend_at_bottom=True,
#                                    legend_at_bottom_columns=2,
#                                    x_labels_major_every=5,
#                                    show_minor_x_labels=False,
#                                    truncate_label=-1,
#                                    style=custom_style)

#     print(outbreak_end)
#     line_chart.title = 'Simulation results - Compartment model'
#     line_chart.x_labels = [str(x) for x in range(outbreak_end)]
#     line_chart.add('latent', results['latent'])
#     line_chart.add('infectious', results['infectious'])
#     line_chart.add('susceptible', [{'value': x, 'color': 'rgba(51, 153, 255, 100)'}
#                                    for x in results['susceptible']])
#     # line_chart.add('infectious_symptomatic_travel', results['symptomatic_travel'],
#                    # color='rgba(230, 0, 0, 70)')
#     # line_chart.add('infectious_symptomatic_no_travel', results['symptomatic_no_travel'])
#     line_chart.add('recovered', results['recovered'])
#     # line_chart.add('infectious_asymptomatic', results['asymptomatic'])
#     line_chart.render_in_browser()
#     line_chart.render_to_file(svg_filepath)

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
        inf = [x[1]['infectious_a'] + x[1]['infectious_t'] + x[1]['infectious_nt'] for x in s]
        rec = [x[1]['recovered'] for x in s]

        outbreak_end = rec.index(max(rec))
        for x in [timesteps, sus, lat, inf, rec]:
            x = x[:outbreak_end+1]

        return {'sus': sus, 'lat': lat, 'inf': inf, 'rec': rec, 'timesteps': timesteps}

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
        writer.writerow(['node_name', 'timestep', 'susceptible', 'latent', 'infectious', 'recovered'])
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


# def test_folium():
#     geo_data = {"type": "FeatureCollection",
#                 "features": []}
#     m = folium.Map([0,3], zoom_start=2)
#     tgj = plugins.TimestampedGeoJson(geo_data)
#     m.add_children(tgj)
#     m.save('folium_test.html')


def compute_commuting_flow(input_file, output_file):
    """
    Commuter flow estimation done using the radiation model of traffic flow.
    See article:
    Simini, Filippo, Marta C. González, Amos Maritan, and Albert-László Barabási. 2012.
    “A Universal Model for Mobility and Migration Patterns.”
    Nature 484 (7392): 96–100.
    """


    g = nx.read_graphml(input_file)

    commuter_ratio = 0.11

    counter = 0
    for i in g.nodes_iter():
        print('computing for node {}'.format(counter))
        counter += 1
        pop_i = g.node[i]['pop']
        neighbors = set(nx.neighbors(g, i))
        for j in neighbors:
            pop_j = g.node[j]['pop']
            other_neighbors = neighbors - set(j)
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
                                             )
    print('writing to file')
    nx.write_graphml(g, output_file)


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


if __name__ == '__main__':
    # G = nx.read_graphml('/data/influenza/rwanda/rwa_net.graphml')
    # prune_edges_with_max_distance(G, 50000)
    # nx.write_graphml(G, '/data/influenza/rwanda/rwa_net_pruned.graphml')
    # input_file = '/home/hugo/data/pygleam/2016-07-15_rwa-net.graphml'
    # output_file = '/home/hugo/data/pygleam/rwa-net_cr.graphml'
    # compute_commuting_flow(input_file, output_file)

    # input_file = '/home/hugo/data/pygleam/rwa-net_cr.graphml'
    # output_file = '/home/hugo/data/pygleam/rwa-net_cr_pruned.graphml'
    # prune_edges_with_min_cr(input_file, output_file, 0.001)

    # input_file = '/home/hugo/Projects/gleam/data/rwa_net.graphml'
    # output_file = '/home/hugo/Projects/gleam/data/base_nodes.jsonp'
    # generate_geojson_base_nodes(input_file, output_file)

    # plot_epidemic_curve_from_results('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
    #                                  'charts/H1N1_node-890_seed-1_n-100_vacc-0.0.svg')

    # get_csv_data_from_results_global('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
    #                                  'csv/H1N1_node-890_seed-1_n-100_vacc-0.0_GLOBAL.csv')

    get_csv_data_from_results_by_node('output/H1N1_node-890_seed-1_n-100_vacc-0.0.jsonp',
                                      'csv/H1N1_node-890_seed-1_n-100_vacc-0.0_BYNODE.csv')