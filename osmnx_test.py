import osmnx as ox
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi,voronoi_plot_2d
import itertools as it
from bentley_ottmann import ensure_planar_graph
from dual_graph import build_dual_graph
import poly_point_isect

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

def generate_coordinates(G, is_multigraph=False, layout_func=None):
    if layout_func is None:
        # Default to spring layout for better edge separation
        layout_func = nx.spring_layout

    G_with_coords = G.copy()
    # Generate positions
    try:
        # Handle different graph types for layout generation
        layout_graph = G.to_undirected() if hasattr(G, 'to_undirected') else G
        if is_multigraph:
            # Convert multigraph to simple graph for layout
            simple_graph = nx.Graph()
            simple_graph.add_nodes_from(layout_graph.nodes())
            for u, v in layout_graph.edges():
                if not simple_graph.has_edge(u, v):
                    simple_graph.add_edge(u, v)
            positions = layout_func(simple_graph)
        else:
            positions = layout_func(layout_graph)
        
        # Add coordinates to nodes
        for node, (x, y) in positions.items():
            G_with_coords.nodes[node]['x'] = x
            G_with_coords.nodes[node]['y'] = y
            print(x, y)

            # nx.draw(G_with_coords, positions, node_size=100)
            
            print(f"Added coordinates using {layout_func.__name__}")
            
    except Exception as e:
        print(f"Error generating layout: {e}")
        # Fallback to simple grid layout
        nodes = list(G.nodes())
        import math
        grid_size = int(math.ceil(math.sqrt(len(nodes))))
        
        for i, node in enumerate(nodes):
            x = i % grid_size
            y = i // grid_size
            G_with_coords.nodes[node]['x'] = float(x)
            G_with_coords.nodes[node]['y'] = float(y)
        
        print("Used fallback grid layout")

    return G_with_coords

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

G_with_coords = generate_coordinates(G, is_multigraph=False, layout_func=nx.spring_layout)
for n in G_with_coords:
    print(n)

for node in G_with_coords.nodes(data=True):
    print(node[1]['x'], node[1]['y'])

for u, v, data in G_with_coords.edges(data=True):
    print(u, v, data)

print(G_with_coords.nodes[0])
nx.draw(G_with_coords, {n: (data['x'], data['y']) for n, data in G_with_coords.nodes(data=True)}, node_size=100)

if not is_planar:
    # for i in range(0, 10):
    #     G_planar, intersections = ensure_planar_graph(G, remove_intersections=False, is_multigraph=False, add_coordinates=True, layout_func=nx.spring_layout)
    #     is_planar, embedding = nx.check_planarity(G_planar.to_undirected())
    #     print(is_planar)
    #     print(embedding)
    coordinates = []
    counter = 0
    for i in range(0, 10):
        node = G_with_coords.nodes(data=True)[counter]
        coordinates.append((float(node['x']), float(node['y'])))
        if i >= 5:
            counter += 2
        else:
            counter += 1
        counter = counter % len(G_with_coords.nodes())
    # poly = [((data['x'], data['y']), (data2['x'], data2['y'])) for u, v, data in G_with_coords.edges(data=True) for data2 in [G_with_coords.nodes[v]]]
    poly: tuple[tuple[float, float], ...] = tuple(tuple(coord) for coord in coordinates)
    print(poly)
    print(type(poly))
    print(type(poly[0]))
    print(type(poly[0][0]))
    intersection_points = poly_point_isect.isect_polygon(poly, validate=True)
    print(intersection_points)
    intersection_point_graph = nx.Graph()
    for i, (x, y) in enumerate(intersection_points):
        intersection_point_graph.add_node(i, pos=(x, y), x=x, y=y)
    nx.draw(intersection_point_graph, {i: (x, y) for i, (x, y) in enumerate(intersection_points)}, node_size=5)
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




