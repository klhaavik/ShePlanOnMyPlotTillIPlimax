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
            # print(x, y)

            # nx.draw(G_with_coords, positions, node_size=100)
            
            # print(f"Added coordinates using {layout_func.__name__}")
            
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

def remove_intersecting_edges(G_planar, intersections, is_multigraph=False):
    removed_edges = []
    edges_to_remove = set()

    X = 0
    Y = 1
    Start = 0
    End = 1

    for intersection_point, (seg1, seg2) in intersections:
        # print(seg1, seg2)
        # Calculate edge lengths
        len1 = ((seg1[End][X] - seg1[Start][X])**2 + (seg1[End][Y] - seg1[Start][Y])**2)**0.5
        len2 = ((seg2[End][X] - seg2[Start][X])**2 + (seg2[End][Y] - seg2[Start][Y])**2)**0.5
        
        # Remove the longer edge
        edge_to_remove = seg2 if len1 < len2 else seg1
        # print(edge_to_remove)
        for node in G_planar.nodes(data=True):
            # print("Graph node coords:", node[1]['x'], node[1]['y'])
            # print("Edge start coords:", edge_to_remove[Start][X], edge_to_remove[Start][Y])
            # print("Edge end coords:", edge_to_remove[End][X], edge_to_remove[End][Y])
            if math.isclose(node[1]['x'], edge_to_remove[Start][X]) and math.isclose(node[1]['y'], edge_to_remove[Start][Y]):
                u = node[0]
                # print("Added as u")
            if math.isclose(node[1]['x'], edge_to_remove[End][X]) and math.isclose(node[1]['y'], edge_to_remove[End][Y]):
                v = node[0]
                # print("Added as v")
            # print("\n")
        if is_multigraph:
            edges_to_remove.add((u, v, 0)) # placeholder key
        else:
            edges_to_remove.add((u, v))
        # print(f"Marked edge {edge_to_remove} for removal due to intersection at {intersection_point}")
    
    # Remove edges from graph based on graph type
    for edge_key in edges_to_remove:
        # print(edge_key)
        if is_multigraph:
            # MultiGraph: edge_key is (u, v, key)
            u, v, key = edge_key
            if G_planar.has_edge(u, v, key):
                G_planar.remove_edge(u, v, key)
                removed_edges.append(edge_key)
        else:
            # Regular Graph: edge_key is (u, v)
            u, v = edge_key
            if G_planar.has_edge(u, v):
                G_planar.remove_edge(u, v)
                removed_edges.append(edge_key)
    
    # print(f"Removed {len(removed_edges)} intersecting edges")
    return G_planar, removed_edges

G = ox.graph.graph_from_point((37.79, -122.407), dist=500, network_type="drive", simplify=True)
# print(type(G))
# node_coords = [[data['y'], data['x']] for node, data in G.nodes(data=True)]
# dual = voronoi_to_networkx(node_coords)
# pos = dict(zip(dual.nodes(), dual.nodes()))
# nx.draw(dual, pos,node_size=5)

# G = nx.complete_graph(5)
# print(type(G))

# G = nx.cycle_graph(3)
# print(type(G))

is_planar, embedding = nx.check_planarity(G.to_undirected())
print(is_planar)
print(embedding)

# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
# ax1.set_visible(True)
# ax2.set_visible(True)
# for ax in (ax1, ax2):
#     limits=plt.axis('on') # turns on axis
#     ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
#     ax.xaxis.set_label_text("X")
#     ax.yaxis.set_label_text("Y")
#     ax.set_aspect('equal', adjustable='box')

# G_with_coords = generate_coordinates(G, is_multigraph=False, layout_func=nx.spring_layout)
G_with_coords = G.copy()
# for n in G_with_coords:
#     print(n)

# counter = 0
# for node in G_with_coords.nodes(data=True):
#     print(node)


# print(G_with_coords.edges(data=True))

# print(G_with_coords.nodes[0])

# nx.draw(G_with_coords, {n: (data['x'], data['y']) for n, data in G_with_coords.nodes(data=True)}, ax=ax1, node_size=100)

if not is_planar:
    # for i in range(0, 10):
    #     G_planar, intersections = ensure_planar_graph(G, remove_intersections=False, is_multigraph=False, add_coordinates=True, layout_func=nx.spring_layout)
    #     is_planar, embedding = nx.check_planarity(G_planar.to_undirected())
    #     print(is_planar)
    #     print(embedding)

    segments_to_remove = []
    for u, v, data in G_with_coords.edges(data=True):
        if 'tunnel' in data and data['tunnel'] == 'yes':
            segments_to_remove.append((u, v))
        if 'bridge' in data and data['bridge'] == 'yes':
            segments_to_remove.append((u, v))

    G_with_coords.remove_edges_from(segments_to_remove)

    fig, ax = ox.plot.plot_graph(
        G_with_coords, bgcolor="k", node_color="blue", node_size=50, edge_linewidth=2, edge_color="#333333"
    )
    # algorithm to add nodes traveled in order for K(5) for input to isect_polygon()
    # coordinates = []
    # counter = 0
    # for i in range(0, 10):
    #     node = G_with_coords.nodes(data=True)[counter]
    #     coordinates.append((float(node['x']), float(node['y'])))
    #     if i >= 5:
    #         counter += 2
    #     else:
    #         counter += 1
    #     counter = counter % len(G_with_coords.nodes())
    segments = []
    for u, v, data in G_with_coords.edges(data=True):
        start_point = (G_with_coords.nodes(data=True)[u]['x'], G_with_coords.nodes(data=True)[u]['y'])
        end_point = (G_with_coords.nodes(data=True)[v]['x'], G_with_coords.nodes(data=True)[v]['y'])
        # print(tuple((start_point, end_point)))
        segments.append(tuple((start_point, end_point)))
    # poly = [((data['x'], data['y']), (data2['x'], data2['y'])) for u, v, data in G_with_coords.edges(data=True) for data2 in [G_with_coords.nodes[v]]]
    # poly = tuple((s, c) for (s, c) in poly)
    # print(segments)
    intersection_points = poly_point_isect.isect_segments_include_segments(segments, validate=True)
    # for pair in intersection_points:
    #     print(pair)



    G_planar, removed_edges = remove_intersecting_edges(G_with_coords, intersection_points, is_multigraph=False)
    is_planar, embedding = nx.check_planarity(G_planar.to_undirected())
    
    # intersection_point_graph = nx.Graph()
    # for i, (x, y) in enumerate(intersection_points):
    #     intersection_point_graph.add_node(i, pos=(x, y), x=x, y=y)
    # nx.draw(intersection_point_graph, {i: (x, y) for i, (x, y) in enumerate(intersection_points)}, ax=ax2, node_size=100)
    # fig, ax = ox.plot.plot_graph(
    #     G_planar, bgcolor="k", node_color="blue", node_size=50, edge_linewidth=2, edge_color="#333333"
    # )
    # nx.draw(G_planar, {n: (data['x'], data['y']) for n, data in G_planar.nodes(data=True)}, ax=ax2, node_size=100)


else:
    G_planar = G


fig, ax = ox.plot.plot_graph(
    G_planar, bgcolor="k", node_color="blue", node_size=50, edge_linewidth=2, edge_color="#333333"
)

dual = build_dual_graph(G_planar, use_coordinates=True, weight_edges=False, return_multidigraph=True)
is_planar, embedding = nx.check_planarity(dual)
print(is_planar)
print(embedding)

max_degree_node = max(dual.degree(), key=lambda x: x[1])[0]
print(f"Removing node {max_degree_node} with degree {dual.degree(max_degree_node)}")
dual.remove_node(max_degree_node)
is_planar, embedding = nx.check_planarity(dual)
print(is_planar)
print(embedding)
    
fig, ax = ox.plot.plot_graph(
    dual, bgcolor="k", node_color="red", node_size=50, edge_linewidth=2, edge_color="#333333"
)

pos = nx.planar_layout(dual)
nx.draw(dual, pos=pos, with_labels=True)
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




