# !/usr/bin/env python
#  encoding: utf-8


import pandas as pd
import matplotlib.pyplot as plt

results = pd.read_json('results.json')

results['step'] = range(1000)

results['infectious'] = results[['symptomatic_no_travel',
                                 'symptomatic_travel',
                                 'asymptomatic']].sum(axis=1)

outbreak_peak = results['infectious'].idxmax()

outbreak_end = results['infectious'][outbreak_peak:].idxmin()
print(outbreak_end)

peak_inf_percentage = results['infectious'].max() / 2000.0
textstr = 'most infectious day: {}\nhighest infection %%: {}%%\nend of outbreak: {}'.format(
            outbreak_peak, peak_inf_percentage, outbreak_end)

fig, ax = plt.subplots(1)

ax.plot(results['recovered'], 'g', label='recovered')
ax.plot(results['susceptible'], 'b', label='susceptible')
ax.plot(results['latent'], 'y', label='latent')
ax.plot(results['infectious'], 'r', label='infectious')
ax.plot(results['asymptomatic'], '0.3', label='asymptomatic')
ax.plot(results['symptomatic_travel'], '0.8', label='symptomatic_travel')
ax.plot(results['symptomatic_no_travel'], '0.5', label='symptomatic_no_travel')
ax.axvline(outbreak_peak, color='r', alpha=0.5, linestyle='--')
ax.axvline(outbreak_end, color='g', alpha=0.5, linestyle='--')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
         verticalalignment='top', bbox=props)
plt.legend()
plt.xlim(0, outbreak_end + 10)


plt.show()
