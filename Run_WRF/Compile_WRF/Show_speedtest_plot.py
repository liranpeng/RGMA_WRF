import matplotlib.pyplot as plt
import numpy as np

# Data from the table
experiments = np.arange(1, 15)
nodes = np.array([4,4,4,4,8,8,8,8,12,12,12,16,16,16])
omp_threads = np.array([1,2,4,8,1,2,4,8,2,4,8,2,4,8])
node_hours = np.array([41.2,32.0,28.8,26.2,72.0,48.0,48.0,48.0,82.3,72.0,67.6,96.0,96.0,87.2])
sim_hours = np.array([10.3,8.0,7.2,6.55,9.0,6.0,6.0,6.0,6.86,6.0,5.63,6.0,6.0,5.45])

# Colors for OMP threads
colors = {1: 'blue', 2: 'green', 4: 'orange', 8: 'red'}
labels = {1: 'OMP=1', 2: 'OMP=2', 4: 'OMP=4', 8: 'OMP=8'}

fig, axs = plt.subplots(1, 2, figsize=(14,6), sharey=False)

# Left plot: Node-hours vs Nodes
for omp in np.unique(omp_threads):
    idx = omp_threads == omp
    # Sort by nodes to connect with lines
    sorted_idx = np.argsort(nodes[idx])
    x = nodes[idx][sorted_idx]
    y = node_hours[idx][sorted_idx]
    axs[0].plot(x, y, marker='o', label=labels[omp], color=colors[omp])
axs[0].set_xlabel("Number of Nodes")
axs[0].set_ylabel("Node-Hours to Simulate 1 Day")
axs[0].set_title("Node-Hours vs Nodes")
axs[0].grid(True)
axs[0].legend(title="OMP Threads")

# Right plot: Sim Hours vs Nodes
for omp in np.unique(omp_threads):
    idx = omp_threads == omp
    sorted_idx = np.argsort(nodes[idx])
    x = nodes[idx][sorted_idx]
    y = sim_hours[idx][sorted_idx]
    axs[1].plot(x, y, marker='o', label=labels[omp], color=colors[omp])
axs[1].set_xlabel("Number of Nodes")
axs[1].set_ylabel("Hours to Simulate 1 Day")
axs[1].set_title("Simulation Hours vs Nodes")
axs[1].grid(True)
axs[1].legend(title="OMP Threads")

plt.tight_layout()
plt.show()

