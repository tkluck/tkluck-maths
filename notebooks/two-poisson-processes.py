import numpy
from scipy import stats

N = 100000

# %% model 1: changes and incidents are independent. Say, a change happens
# on average every two units of time; an incident happens every three units
# of time
changes = stats.expon(scale=2).rvs(N)
incidents = stats.expon(scale=3).rvs(N)

change_times   = numpy.cumsum(changes)
incident_times = numpy.cumsum(incidents)

# %% Restrict to the first few instances
change_times = [t for t in change_times if t < N/10]
incident_times = [t for t in incident_times if t < N/10]
timeline = [(t, "change") for t in change_times]
timeline += [(t, "incident") for t in incident_times]
timeline.sort()

# %% scatter plot to see what's happening
import matplotlib.pyplot as plt
plt.scatter(change_times, [0 for _ in change_times], color="lightblue")
plt.scatter(incident_times, [0 for _ in incident_times], color="orange")
plt.xlim(xmin=0, xmax=40)

# %%
deltas = []
for ((t1, w1), (t2,w2)) in zip(timeline[:-1], timeline[1:]):
    if w2 == "incident" and w1 == "change":
        deltas.append(t2-t1)

plt.hist(deltas, bins=50);

# %% model 2: two-thirds of changes are followed by an incident, about 0.3
# units of time later. Note that this means that an incident occurs, on average
# every 3 units of time, just like in model 1.
# (run this cell, and then run all the cells above again to see the change in the plot.)

changes        = stats.expon(scale=2).rvs(N)
incident_delays= stats.expon(scale=0.3).rvs(N)
change_times   = numpy.cumsum(changes)
has_incident   = 3*stats.uniform().rvs(N) < 2

incident_times = (change_times + incident_delays)[has_incident]

# %% validate the claim about average number of incidents

numpy.average(incident_times[1:] - incident_times[:-1])
