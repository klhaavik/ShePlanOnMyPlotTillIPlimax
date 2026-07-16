import networkx as nx
from collections import defaultdict

def get_faces(G):
    """
    Extract all faces from a planar graph.
    
    Args:
        G: A planar NetworkX graph
        
    Returns:
        List of faces, where each face is a list of vertices in order around the face
    """
    if not G.edges():
        # Handle trivial cases
        if not G.nodes():
            return []
        else:
            return [list(G.nodes())]  # Single face containing all isolated vertices
    
    # Get planar embedding - this gives us the cyclic order of neighbors
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        raise ValueError("Graph is not planar")
    
    faces = []
    used_half_edges = set()
    
    # For each edge, we have two "half-edges" (u,v) and (v,u)
    # Each half-edge borders exactly one face
    for u in G.nodes():
        # Get neighbors in the planar embedding order
        neighbors = list(embedding[u])
        
        for i, v in enumerate(neighbors):
            half_edge = (u, v)
            
            if half_edge in used_half_edges:
                continue
                
            # Trace the face starting from this half-edge
            face = []
            current = u
            next_node = v
            
            while True:
                face.append(current)
                used_half_edges.add((current, next_node))
                
                # Move to the next node
                current = next_node
                
                # Find the next edge in the face by going clockwise around current
                current_neighbors = list(embedding[current])
                
                # Find where we came from
                prev_idx = current_neighbors.index(face[-1])
                
                # Next edge is the one after the previous edge (clockwise)
                next_idx = (prev_idx + 1) % len(current_neighbors)
                next_node = current_neighbors[next_idx]
                
                # If we've completed the cycle, break
                if current == u and next_node == v:
                    break
                    
                # Safety check to prevent infinite loops
                if len(face) > len(G.edges()) + len(G.nodes()):
                    break
            
            if len(face) >= 3:  # Valid face must have at least 3 vertices
                faces.append(face)
    
    return faces


def calculate_face_centroid(face, pos):
    """Calculate the centroid of a face given vertex positions"""
    if not face or not pos:
        return (0, 0)
    
    x_coords = [pos[v][0] for v in face if v in pos]
    y_coords = [pos[v][1] for v in face if v in pos]
    
    if not x_coords or not y_coords:
        return (0, 0)
    
    centroid_x = sum(x_coords) / len(x_coords)
    centroid_y = sum(y_coords) / len(y_coords)
    
    return (centroid_x, centroid_y)


def calculate_edge_midpoint(u, v, pos):
    """Calculate the midpoint of an edge given vertex positions"""
    if u not in pos or v not in pos:
        return (0, 0)
    
    x1, y1 = pos[u]
    x2, y2 = pos[v]
    
    midpoint_x = (x1 + x2) / 2
    midpoint_y = (y1 + y2) / 2
    
    return (midpoint_x, midpoint_y)


def euclidean_distance(p1, p2):
    """Calculate Euclidean distance between two points"""
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5


def build_dual_graph(G, use_coordinates=True, weight_edges=True, return_multidigraph=False):
    """
    Build the dual graph of a planar graph G with coordinate system support.
    
    Args:
        G: A planar NetworkX graph
        use_coordinates: If True, use node coordinates for dual node positioning
        weight_edges: If True, weight dual edges by distance between face centroids
        return_multidigraph: If True, return MultiDiGraph with directed edges for each original edge
        
    Returns:
        NetworkX graph representing the dual graph (Graph, DiGraph, or MultiDiGraph)
    """
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        raise ValueError("Graph is not planar")

    faces = get_faces(G)
    # print(faces)
    
    # Choose the appropriate graph type
    if return_multidigraph:
        dual = nx.MultiDiGraph()
    else:
        dual = nx.Graph()
    
    # Get node positions if they exist
    pos = None
    if use_coordinates and hasattr(G, 'nodes'):
        # Try to get positions from various common attributes
        pos_attrs = ['pos', 'position', 'coordinates', 'coord']
        for attr in pos_attrs:
            if all(attr in G.nodes[node] for node in G.nodes() if G.nodes[node]):
                pos = {node: G.nodes[node][attr] for node in G.nodes()}
                break
        
        # Also check if nodes have x,y attributes
        if pos is None and all('x' in G.nodes[node] and 'y' in G.nodes[node] 
                              for node in G.nodes() if G.nodes[node]):
            pos = {node: (G.nodes[node]['x'], G.nodes[node]['y']) for node in G.nodes()}

    # Each face corresponds to a node in the dual
    dual_pos = {}
    for i, face in enumerate(faces):
        # Store the original face information
        dual.add_node(i, face=face, face_size=len(face))
        
        # Calculate face centroid if coordinates are available
        if pos and use_coordinates:
            x, y = calculate_face_centroid(face, pos)
            dual.nodes[i]['x'] = x
            dual.nodes[i]['y'] = y
            dual_pos[i] = (x, y)
            
            # Calculate face area (for polygons)
            if len(face) >= 3:
                area = calculate_polygon_area(face, pos)
                dual.nodes[i]['area'] = area

    # Set dual graph position attribute if coordinates were used
    if dual_pos and use_coordinates:
        nx.set_node_attributes(dual, dual_pos, 'pos')

    if return_multidigraph:
        # For multidigraph: each original edge becomes a directed edge in the dual
        # Direction is determined by the orientation of the edge relative to each face
        
        for u, v in G.edges():
            edge = frozenset((u, v))
            faces_containing_edge = []
            
            # Find which faces contain this edge and in what orientation
            for i, face in enumerate(faces):
                n = len(face)
                for j in range(n):
                    face_u = face[j]
                    face_v = face[(j + 1) % n]
                    
                    if (face_u == u and face_v == v) or (face_u == v and face_v == u):
                        # Determine the orientation relative to this face
                        if face_u == u and face_v == v:
                            # Edge goes u->v in this face (counterclockwise orientation)
                            faces_containing_edge.append((i, 'ccw'))
                        else:
                            # Edge goes v->u in this face (clockwise orientation)
                            faces_containing_edge.append((i, 'cw'))
                        break
            
            # Add directed edges in the dual
            if len(faces_containing_edge) == 2:
                (f1, orient1), (f2, orient2) = faces_containing_edge
                
                # Create edge attributes
                edge_attrs = {
                    'original_edge': (u, v),
                    'original_edge_set': edge
                }
                
                # Calculate edge weight based on distance between face centroids
                if pos and use_coordinates and weight_edges and f1 in dual_pos and f2 in dual_pos:
                    distance = euclidean_distance(dual_pos[f1], dual_pos[f2])
                    edge_attrs['weight'] = distance
                    edge_attrs['distance'] = distance
                    
                    # Also store the midpoint of the original edge
                    midpoint = calculate_edge_midpoint(u, v, pos)
                    edge_attrs['edge_midpoint'] = midpoint
                
                # Add two directed edges (one for each direction of the original edge)
                # The direction in the dual is perpendicular to the original edge
                if orient1 == 'ccw':
                    # If edge is ccw in face f1, then dual edge goes from f1 to f2
                    dual.add_edge(f1, f2, **{**edge_attrs, 'orientation': f'{u}->{v}'})
                else:
                    # If edge is cw in face f1, then dual edge goes from f2 to f1
                    dual.add_edge(f2, f1, **{**edge_attrs, 'orientation': f'{u}->{v}'})
                
    else:
        # Original behavior for simple graphs
        edge_faces = {}
        
        # Map each edge to the faces it bounds
        for i, face in enumerate(faces):
            n = len(face)
            for j in range(n):
                u = face[j]
                v = face[(j + 1) % n]
                edge = frozenset((u, v))
                edge_faces.setdefault(edge, []).append(i)

        # For each edge in original graph, connect the faces that share it
        for edge, face_indices in edge_faces.items():
            if len(face_indices) == 2:
                f1, f2 = face_indices
                
                # Add edge attributes
                edge_attrs = {'shared_edge': edge}
                
                # Calculate edge weight based on distance between face centroids
                if pos and use_coordinates and weight_edges and f1 in dual_pos and f2 in dual_pos:
                    distance = euclidean_distance(dual_pos[f1], dual_pos[f2])
                    edge_attrs['weight'] = distance
                    edge_attrs['distance'] = distance
                    
                    # Also store the midpoint of the original edge
                    u, v = edge
                    midpoint = calculate_edge_midpoint(u, v, pos)
                    edge_attrs['edge_midpoint'] = midpoint
                
                dual.add_edge(f1, f2, **edge_attrs)

    # Add metadata about the coordinate system if available
    if pos and use_coordinates:
        dual.graph['has_coordinates'] = True
        dual.graph['crs'] = G.graph["crs"]  # Could be extended for other systems
        
        # Calculate bounding box of original graph
        x_coords = [pos[v][0] for v in pos]
        y_coords = [pos[v][1] for v in pos]
        dual.graph['bbox'] = {
            'min_x': min(x_coords),
            'max_x': max(x_coords),
            'min_y': min(y_coords),
            'max_y': max(y_coords)
        }
    else:
        dual.graph['has_coordinates'] = False
    
    # Add metadata about the dual type
    dual.graph['dual_type'] = 'multidigraph' if return_multidigraph else 'simple'
    dual.graph['original_edges'] = len(G.edges())
    dual.graph['dual_edges'] = len(dual.edges())

    return dual


def calculate_polygon_area(vertices, pos):
    """Calculate the area of a polygon using the shoelace formula"""
    if len(vertices) < 3:
        return 0.0
    
    n = len(vertices)
    area = 0.0
    
    for i in range(n):
        j = (i + 1) % n
        if vertices[i] in pos and vertices[j] in pos:
            x1, y1 = pos[vertices[i]]
            x2, y2 = pos[vertices[j]]
            area += x1 * y2 - x2 * y1
    
    return abs(area) / 2.0


# Example usage and test
if __name__ == "__main__":
    # Test with a simple triangle with coordinates
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0)])
    
    # Add coordinate positions
    pos = {0: (0, 0), 1: (1, 0), 2: (0.5, 0.866)}  # Equilateral triangle
    nx.set_node_attributes(G, pos, 'pos')
    
    faces = get_faces(G)
    print("Faces of triangle:", faces)
    
    # Test simple dual
    dual_simple = build_dual_graph(G, use_coordinates=True, weight_edges=True, return_multidigraph=False)
    print("\nSimple dual graph nodes:", dual_simple.nodes(data=True))
    print("Simple dual graph edges:", dual_simple.edges(data=True))
    
    # Test multidigraph dual
    dual_multi = build_dual_graph(G, use_coordinates=True, weight_edges=True, return_multidigraph=True)
    print("\nMultidigraph dual nodes:", dual_multi.nodes(data=True))
    print("Multidigraph dual edges:", list(dual_multi.edges(data=True)))
    print("Multidigraph metadata:", dual_multi.graph)
    
    # Test with a square (more interesting for multidigraph)
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])  # Square with diagonal
    
    # Add coordinate positions for a unit square
    pos2 = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1)}
    nx.set_node_attributes(G2, pos2, 'pos')
    
    print(f"\nOriginal graph G2 has {len(G2.edges())} edges")
    
    dual2_multi = build_dual_graph(G2, use_coordinates=True, weight_edges=True, return_multidigraph=True)
    print(f"Multidigraph dual has {len(dual2_multi.edges())} directed edges")
    print("Square with diagonal - multidigraph dual edges:")
    for u, v, key, data in dual2_multi.edges(keys=True, data=True):
        print(f"  {u} -> {v} (key={key}): {data.get('original_edge', '')}, orientation: {data.get('orientation', '')}")
    
    # Demonstrate the difference
    dual2_simple = build_dual_graph(G2, use_coordinates=True, weight_edges=True, return_multidigraph=False)
    print(f"Simple dual has {len(dual2_simple.edges())} undirected edges")
    
    print(f"\nEdge count comparison for square with diagonal:")
    print(f"  Original graph: {len(G2.edges())} edges")
    print(f"  Simple dual: {len(dual2_simple.edges())} edges") 
    print(f"  Multidigraph dual: {len(dual2_multi.edges())} directed edges")