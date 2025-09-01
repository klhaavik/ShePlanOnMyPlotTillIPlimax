import osmnx as ox
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi,voronoi_plot_2d
import itertools as it
from bentley_ottmann import ensure_planar_graph
from dual_graph import build_dual_graph

def voronoi_to_networkx(points):
# we get the voronoi diagram
  vor = Voronoi(points)
  voronoi_plot_2d(vor)
  print(vor.ridge_vertices)

  G = nx.Graph()

# Add an edge for each ridge in the Voronoi diagram that connects two points in the range [0,1] 
  for simplex in vor.ridge_vertices:
      if -1 not in simplex:
          i, j = simplex
          p = vor.vertices[i]
          q = vor.vertices[j]
          if 0 <= p[0] <= 1 and 0 <= p[1] <= 1 and 0 <= q[0] <= 1 and 0 <= q[1] <= 1:
              distance = np.linalg.norm(p - q) # Calculate the Euclidean distance between p and q
              G.add_edge(tuple(p), tuple(q),weight=distance)

  return G

def draw_graphs_with_coords(graphs):
    # get geopandas dataframes
    gdfs = []
    for i in range(len(graphs)):
        gdfs[i] = ox.graph_to_gdfs(graphs[i]) # returns a tuple [gdf_nodes, gdf_edges]

    # create a figure with two plots next to each other
    fig, axs = plt.subplots(nrows=1, ncols=len(graphs))

    for i in range(len(graphs)):
        gdfs[i].plot(ax=axs[i], linewidth=0.1)

def draw_graphs(graphs, node_size=5):
    fig, axs = plt.subplots(nrows=1, ncols=len(graphs))

    for i in range(len(graphs)):
        pos = nx.spring_layout(graphs[i])
        nx.draw(graphs[i], pos, node_size=node_size, ax=axs[i])

G = ox.graph.graph_from_point((37.79, -122.41), dist=750, network_type="drive", simplify=True)
print(type(G))
# node_coords = [[data['y'], data['x']] for node, data in G.nodes(data=True)]
# dual = voronoi_to_networkx(node_coords)
# pos = dict(zip(dual.nodes(), dual.nodes()))
# nx.draw(dual, pos,node_size=5)

G = nx.complete_graph(5)
print(type(G))

# G = nx.cycle_graph(3)
# print(type(G))

is_planar, embedding = nx.check_planarity(G.to_undirected())
print(is_planar)
print(embedding)

if not is_planar:
    G_planar, intersections = ensure_planar_graph(G, remove_intersections=False, is_multigraph=False, add_coordinates=True, layout_func=nx.circular_layout)
    is_planar, embedding = nx.check_planarity(G_planar.to_undirected())
    print(is_planar)
    print(embedding)
else:
    G_planar = G


# fig, ax = ox.plot.plot_graph(
#     G_planar, bgcolor="k", node_color="blue", node_size=50, edge_linewidth=2, edge_color="#333333"
# )

# dual = build_dual_graph(G_planar, use_coordinates=False, return_multidigraph=True)
# is_planar, embedding = nx.check_planarity(dual)
# print(is_planar)
# print(embedding)
# print(dual.nodes(data=True))

# G_with_coords = G.copy()
# positions = nx.spring_layout(G)

# for node, (x, y) in positions.items():
#     G_with_coords.nodes[node]['x'] = x
#     G_with_coords.nodes[node]['y'] = y
#     print(x, y)

# nx.draw(G_with_coords, positions, node_size=100)

# draw_graphs_with_coords([G_planar, dual])
# draw_graphs([G_planar, dual], node_size=100)



# fig, ax = ox.plot.plot_graph(
#     dual, ax=bgcolor="k", node_color="blue", node_size=50, edge_linewidth=2, edge_color="#333333"
# )

def draw_labeled_multigraph(G, attr_name, ax=None):
    """
    Length of connectionstyle must be at least that of a maximum number of edges
    between pair of nodes. This number is maximum one-sided connections
    for directed graph and maximum total connections for undirected graph.
    """

    pos = nx.shell_layout(G)
    
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=20, ax=ax)
    
    nx.draw_networkx_edges(G, pos, edge_color="grey", ax=ax)

    labels = {
        (u, v): f"{attr_name}={attrs[attr_name]}" 
        for u, v, attrs in G.edges(data=True)
        if attr_name in attrs
    }

    nx.draw_networkx_edge_labels(
        G,
        pos,
        labels,
        label_pos=0.3,
        font_color="blue",
        bbox={"alpha": 0},
        ax=ax,
    )

# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
# draw_labeled_multigraph(G, "w", ax1)

# dual = build_dual_graph(G)
# is_planar, embedding = nx.check_planarity(dual)
# print(is_planar)
# print(embedding)
# draw_labeled_multigraph(dual, "w", ax2)
plt.show()


def faces_sharing_edge(edge, faces):
    sharing_faces = []
    for face_id, face_edges in faces.items():
        if edge in face_edges or (edge[1], edge[0]) in face_edges:
            sharing_faces.append(face_id)
    return sharing_faces

# ox.io.save_graph_geopackage(G, filepath="./san_fran_default.gpkg")
# dual = nx.MultiGraph()
# for face_id, face_edges in faces.items():
#     dual.add_node(face_id)

# for edge in G.edges():
#     face1, face2 = faces_sharing_edge(edge, faces)
#     dual.add_edge(face1, face2, original_edge=edge)




