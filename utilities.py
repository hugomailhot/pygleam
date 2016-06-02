# !/usr/bin/env python
#  encoding: utf-8

import pandas as pd
import pygal

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

    custom_style = pygal.style.Style(colors=["#fdae61","#d7191c","#a6d96a","#1a9641"])
    custom_style = pygal.style.DarkStyle
    custom_style.colors = ["#feed6c","#bf4646","#56c2d6","#516083"]
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
