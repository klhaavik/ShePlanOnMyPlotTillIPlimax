import osmnx as ox
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi,voronoi_plot_2d
import itertools as it
from bentley_ottmann import ensure_planar_graph
from dual_graph import build_dual_graph
import poly_point_isect
import math

G = nx.Graph()
G.add_nodes_from([1,2,3,4,5,6])
G.add_edges_from([(1,2),(2,3),(3,4),(4,5),(5,1),(2,6)])
nx.draw(G, with_labels=True)
plt.show()

dual = build_dual_graph(G, use_coordinates=False, return_multidigraph=True)
is_planar, embedding = nx.check_planarity(dual)
print(is_planar)
nx.draw(dual, with_labels=True)
plt.show()
