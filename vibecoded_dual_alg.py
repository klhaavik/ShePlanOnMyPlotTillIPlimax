import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Tuple, Any

def edge_to_vertex_dual(G: nx.Graph) -> nx.Graph:
    """
    Convert a street network graph to its edge-to-vertex dual (line graph).
    
    In the dual:
    - Each edge in the original becomes a node
    - Two nodes in the dual are connected if their corresponding edges 
      in the original share a vertex
    - Node positions are midpoints of their corresponding original edges
    
    Args:
        G: NetworkX graph with node attributes 'x' and 'y' for coordinates
        
    Returns:
        NetworkX graph representing the edge-to-vertex dual
    """
    # Create the dual graph
    dual = nx.Graph() if not G.is_directed() else nx.DiGraph()
    
    # Map from original edges to dual nodes
    edge_to_node = {}
    node_counter = 0
    
    # Step 1: Create dual nodes from original edges
    for edge in G.edges():
        u, v = edge
        
        # Get coordinates of the edge endpoints
        u_x, u_y = G.nodes[u]['x'], G.nodes[u]['y']
        v_x, v_y = G.nodes[v]['x'], G.nodes[v]['y']
        
        # Calculate midpoint coordinates for the dual node
        mid_x = (u_x + v_x) / 2
        mid_y = (u_y + v_y) / 2
        
        # Add node to dual graph
        dual_node_id = f"e_{u}_{v}"
        dual.add_node(dual_node_id, x=mid_x, y=mid_y, 
                     original_edge=(u, v))
        
        # Store mapping
        edge_to_node[edge] = dual_node_id
        if not G.is_directed():
            # For undirected graphs, also map the reverse edge
            edge_to_node[(v, u)] = dual_node_id
    
    # Step 2: Create dual edges between nodes whose original edges share a vertex
    original_edges = list(G.edges())
    
    for i, edge1 in enumerate(original_edges):
        for j, edge2 in enumerate(original_edges[i+1:], i+1):
            # Check if edges share a vertex
            shared_vertices = set(edge1) & set(edge2)
            
            if shared_vertices:
                dual_node1 = edge_to_node[edge1]
                dual_node2 = edge_to_node[edge2]
                
                # Add edge in dual graph
                dual.add_edge(dual_node1, dual_node2)
    
    return dual

def visualize_dual_transformation(G: nx.Graph, figsize=(15, 6)):
    """
    Visualize the original graph and its edge-to-vertex dual side by side.
    """
    dual_G = edge_to_vertex_dual(G)
    for u, v, data in dual_G.edges(data=True):
        print(f"Original edge: ({u}, {v}) with data {data}")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Plot original graph
    pos_orig = {node: (G.nodes[node]['x'], G.nodes[node]['y']) for node in G.nodes()}
    nx.draw(G, pos_orig, ax=ax1, with_labels=True, node_color='lightblue', 
            node_size=500, font_size=8, font_weight='bold')
    ax1.set_title("Original Street Network\n(intersections as nodes, streets as edges)")
    ax1.set_aspect('equal')
    
    # Plot dual graph
    pos_dual = {node: (dual_G.nodes[node]['x'], dual_G.nodes[node]['y']) 
                for node in dual_G.nodes()}
    nx.draw(dual_G, pos_dual, ax=ax2, with_labels=True, node_color='lightcoral', 
            node_size=1000, font_size=6, font_weight='bold')
    ax2.set_title("Edge-to-Vertex Dual\n(streets as nodes, connections as edges)")
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()

# Example usage and testing
def create_sample_street_network():
    """Create a sample street network for testing."""
    G = nx.MultiDiGraph()
    
    # Add intersections with coordinates
    intersections = {
        'A': (0, 0),
        'B': (2, 0),
        'C': (4, 0),
        'D': (1, 2),
        'E': (3, 2)
    }
    
    for node, (x, y) in intersections.items():
        G.add_node(node, x=x, y=y)
    
    # Add streets (edges) - some bidirectional (add both directions)
    streets = [
        ('A', 'B'), ('B', 'A'),  # bidirectional street
        ('B', 'C'), ('C', 'B'),  # bidirectional street
        ('A', 'D'),              # one-way street
        ('B', 'D'), ('D', 'B'),  # bidirectional street
        ('B', 'E'),              # one-way street
        ('C', 'E'), ('E', 'C'),  # bidirectional street
        ('D', 'E'), ('E', 'D')   # bidirectional street
    ]
    
    G.add_edges_from(streets)
    return G

# Test the implementation
if __name__ == "__main__":
    # Create sample network
    print("Creating sample street network...")
    street_network = create_sample_street_network()
    
    print(f"Original network: {len(street_network.nodes())} intersections, {len(street_network.edges())} streets")
    
    # Generate dual
    dual_network = edge_to_vertex_dual(street_network)
    print(f"Dual network: {len(dual_network.nodes())} nodes, {len(dual_network.edges())} edges")
    
    # Show some dual node details
    print("\nDual node coordinates (street midpoints):")
    for node in list(dual_network.nodes())[:3]:  # Show first 3
        x, y = dual_network.nodes[node]['x'], dual_network.nodes[node]['y']
        orig_edge = dual_network.nodes[node]['original_edge']
        print(f"  {node}: ({x:.1f}, {y:.1f}) - midpoint of edge {orig_edge}")
    
    # Visualize
    print("\nGenerating visualization...")
    visualize_dual_transformation(street_network)