# !/usr/bin/env python
#  encoding: utf-8

import pandas as pd
import pygal
import networkx as nx
import json
import folium
from folium import plugins

def plot_results(res, svg_filepath):
    results = pd.DataFrame(res)
    results['infectious'] = results[['infectious_nt',
                                     'infectious_t',
                                     'infectious_a']].sum(axis=1)

    outbreak_peak = results['infectious'].idxmax()
    outbreak_end = results['infectious'][outbreak_peak:].idxmin()
    peak_inf_percentage = results['infectious'].max() / 2000.0

    textstr = ('outbreak peak: {}\n'.format(outbreak_peak) +
               'highest infection %%: {}%%\n'.format(peak_inf_percentage) +
               'outbreak end: {}'.format(outbreak_end))

    results = results[:outbreak_end]

    custom_style = pygal.style.Style(colors=["#fdae61", "#d7191c", "#a6d96a", "#1a9641"])
    custom_style = pygal.style.DarkStyle
    custom_style.colors = ["#feed6c", "#bf4646", "#56c2d6", "#516083"]
    line_chart = pygal.StackedLine(dots_size=1,
                                   fill=True,
                                   show_dots=False,
                                   show_y_guides=False,
                                   legend_at_bottom=True,
                                   legend_at_bottom_columns=2,
                                   x_labels_major_every=5,
                                   show_minor_x_labels=False,
                                   truncate_label=-1,
                                   style=custom_style)

    print(outbreak_end)
    line_chart.title = 'Simulation results - Compartment model'
    line_chart.x_labels = [str(x) for x in range(outbreak_end)]
    line_chart.add('latent', results['latent'])
    line_chart.add('infectious', results['infectious'])
    line_chart.add('susceptible', [{'value': x, 'color': 'rgba(51, 153, 255, 100)'}
                                   for x in results['susceptible']])
    # line_chart.add('infectious_symptomatic_travel', results['symptomatic_travel'],
                   # color='rgba(230, 0, 0, 70)')
    # line_chart.add('infectious_symptomatic_no_travel', results['symptomatic_no_travel'])
    line_chart.add('recovered', results['recovered'])
    # line_chart.add('infectious_asymptomatic', results['asymptomatic'])
    line_chart.render_in_browser()
    line_chart.render_to_file(svg_filepath)


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
    graph = nx.read_graphml(input_filepath)
    bad_edges = [x for x in graph.edges_iter(data=True)
                 if x[2]['commuting_rate'] < min_cr]

    graph.remove_edges_from(bad_edges)
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


def test_folium():
    geo_data = {"type": "FeatureCollection",
                "features": []}
    m = folium.Map([0,3], zoom_start=2)
    tgj = plugins.TimestampedGeoJson(geo_data)
    m.add_children(tgj)
    m.save('folium_test.html')


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

    for i in g.nodes_iter():
        pop_i = g.node[i]['pop']
        neighbors = nx.neighbors(g, i)
        for j in neighbors:
            pop_j = g.node[j]['pop']
            other_neighbors = set(neighbors) - set(j)
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


    nx.write_graphml(g, output_file)


def generate_geojson_base_nodes(input_file, output_file):
    """
    Given a network file, will collect latitude and longitude for every node and
    generate a geojson file to allow mapping of the nodes on a Leaflet map.
    """
    def generate_point_object(node):
        return  {
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
    # input_file = '/data/influenza/rwanda/rwa_net_full.graphml'
    # output_file = '/data/influenza/rwanda/rwa_net_full_cr.graphml'
    # compute_commuting_flow(input_file, output_file)
    
    # input_file = '/data/influenza/rwanda/rwa_net_full_cr.graphml'
    # output_file = '/data/influenza/rwanda/rwa_net_full_cr_pruned.graphml'
    # prune_edges_with_min_cr(input_file, output_file, 0.001)

    input_file = '/home/hugo/Projects/gleam/data/rwa_net.graphml'
    output_file = '/home/hugo/Projects/gleam/data/base_nodes.jsonp'
    generate_geojson_base_nodes(input_file, output_file)
