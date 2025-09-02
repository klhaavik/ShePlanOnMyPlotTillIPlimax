import networkx as nx
import osmnx as ox
from collections import namedtuple
from typing import List, Tuple, Set, Optional
import heapq
from enum import Enum
from functools import total_ordering

# Data structures for the algorithm
Point = namedtuple('Point', ['x', 'y'])

class EventType(Enum):
    START = 1
    END = 2
    INTERSECTION = 3

@total_ordering
class Event:
    def __init__(self, x: float, y: float, event_type: EventType, segments: List, event_id: int = 0):
        self.x = x
        self.y = y
        self.event_type = event_type
        self.segments = segments
        self.event_id = event_id  # Unique identifier for tie-breaking
    
    def __lt__(self, other):
        # Primary sort by x-coordinate
        if self.x != other.x:
            return self.x < other.x
        # Secondary sort by y-coordinate
        if self.y != other.y:
            return self.y < other.y
        # Tertiary sort by event type (START < INTERSECTION < END)
        if self.event_type != other.event_type:
            return self.event_type.value < other.event_type.value
        # Final tie-breaker using unique event_id
        return self.event_id < other.event_id
    
    def __eq__(self, other):
        return (self.x == other.x and self.y == other.y and 
                self.event_type == other.event_type and self.event_id == other.event_id)

class Segment:
    def __init__(self, start: Point, end: Point, edge_key: Tuple):
        # Ensure start point has smaller x-coordinate (or smaller y if x is equal)
        if start.x < end.x or (start.x == end.x and start.y < end.y):
            self.start = start
            self.end = end
        else:
            self.start = end
            self.end = start
        self.edge_key = edge_key  # (u, v, key) for networkx multigraph
    
    def __eq__(self, other):
        return self.start == other.start and self.end == other.end
    
    def __hash__(self):
        return hash((self.start, self.end))
    
    def y_at_x(self, x: float) -> float:
        """Calculate y-coordinate at given x on the line segment"""
        if self.start.x == self.end.x:  # Vertical line
            return self.start.y
        
        # Linear interpolation
        t = (x - self.start.x) / (self.end.x - self.start.x)
        return self.start.y + t * (self.end.y - self.start.y)
    
    def intersects(self, other: 'Segment') -> Optional[Point]:
        """Find intersection point between two segments using line intersection formula"""
        x1, y1 = self.start.x, self.start.y
        x2, y2 = self.end.x, self.end.y
        x3, y3 = other.start.x, other.start.y
        x4, y4 = other.end.x, other.end.y
        
        # Calculate denominators
        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        
        if abs(denom) < 1e-10:  # Lines are parallel
            return None
        
        # Calculate intersection parameters
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
        u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
        
        # Check if intersection is within both segments
        if 0 <= t <= 1 and 0 <= u <= 1:
            ix = x1 + t * (x2 - x1)
            iy = y1 + t * (y2 - y1)
            return Point(ix, iy)
        
        return None

class StatusStructure:
    """Red-black tree simulation using sorted list for active segments"""
    def __init__(self):
        self.segments = []
        self.current_x = 0
    
    def update_sweep_line(self, x: float):
        self.current_x = x
        # Re-sort segments based on y-coordinate at current x
        self.segments.sort(key=lambda s: s.y_at_x(x))
    
    def insert(self, segment: Segment):
        self.segments.append(segment)
        self.segments.sort(key=lambda s: s.y_at_x(self.current_x))
    
    def remove(self, segment: Segment):
        if segment in self.segments:
            self.segments.remove(segment)
    
    def get_neighbors(self, segment: Segment) -> Tuple[Optional[Segment], Optional[Segment]]:
        """Get segments immediately above and below the given segment"""
        if segment not in self.segments:
            return None, None
        
        idx = self.segments.index(segment)
        above = self.segments[idx + 1] if idx + 1 < len(self.segments) else None
        below = self.segments[idx - 1] if idx > 0 else None
        return above, below

def bentley_ottmann_intersections(segments: List[Segment]) -> Set[Tuple[Segment, Segment, Point]]:
    """
    Bentley-Ottmann algorithm to find all intersections between line segments
    Returns set of (segment1, segment2, intersection_point) tuples
    """
    if not segments:
        return set()
    
    # Create event queue with unique event IDs
    events = []
    event_counter = 0
    
    # Add start and end events for each segment
    for segment in segments:
        start_event = Event(segment.start.x, segment.start.y, EventType.START, [segment], event_counter)
        heapq.heappush(events, start_event)
        event_counter += 1
        
        end_event = Event(segment.end.x, segment.end.y, EventType.END, [segment], event_counter)
        heapq.heappush(events, end_event)
        event_counter += 1
    
    status = StatusStructure()
    intersections = set()
    processed_intersections = set()
    
    while events:
        event = heapq.heappop(events)
        x, y = event.x, event.y
        event_type = event.event_type
        event_segments = event.segments
        
        status.update_sweep_line(x)
        
        if event_type == EventType.START:
            segment = event_segments[0]
            status.insert(segment)
            
            # Check for intersections with neighbors
            above, below = status.get_neighbors(segment)
            
            if above:
                intersection = segment.intersects(above)
                if intersection and intersection.x > x:
                    key = tuple(sorted([id(segment), id(above)]))
                    if key not in processed_intersections:
                        intersection_event = Event(intersection.x, intersection.y, 
                                                 EventType.INTERSECTION, [segment, above], event_counter)
                        heapq.heappush(events, intersection_event)
                        event_counter += 1
                        processed_intersections.add(key)
            
            if below:
                intersection = segment.intersects(below)
                if intersection and intersection.x > x:
                    key = tuple(sorted([id(segment), id(below)]))
                    if key not in processed_intersections:
                        intersection_event = Event(intersection.x, intersection.y, 
                                                 EventType.INTERSECTION, [segment, below], event_counter)
                        heapq.heappush(events, intersection_event)
                        event_counter += 1
                        processed_intersections.add(key)
        
        elif event_type == EventType.END:
            segment = event_segments[0]
            above, below = status.get_neighbors(segment)
            status.remove(segment)
            
            # Check if newly adjacent segments intersect
            if above and below:
                intersection = above.intersects(below)
                if intersection and intersection.x > x:
                    key = tuple(sorted([id(above), id(below)]))
                    if key not in processed_intersections:
                        intersection_event = Event(intersection.x, intersection.y, 
                                                 EventType.INTERSECTION, [above, below], event_counter)
                        heapq.heappush(events, intersection_event)
                        event_counter += 1
                        processed_intersections.add(key)
        
        elif event_type == EventType.INTERSECTION:
            seg1, seg2 = event_segments
            intersection_point = Point(x, y)
            intersections.add((seg1, seg2, intersection_point))
            
            # Swap segments in status structure
            if seg1 in status.segments and seg2 in status.segments:
                idx1 = status.segments.index(seg1)
                idx2 = status.segments.index(seg2)
                status.segments[idx1], status.segments[idx2] = status.segments[idx2], status.segments[idx1]
    
    return intersections

def ensure_planar_graph(G, remove_intersections: bool = True, is_multigraph: bool = None, 
                       add_coordinates: bool = True, layout_func=None) -> Tuple:
    """
    Ensure graph planarity using Bentley-Ottmann algorithm
    
    Args:
        G: NetworkX graph (typically from OSMnx)
        remove_intersections: If True, remove intersecting edges to make graph planar
        is_multigraph: If None, auto-detect; if True, treat as MultiGraph/MultiDiGraph;
                      if False, treat as Graph/DiGraph
        add_coordinates: If True, add coordinates to graphs that don't have them
        layout_func: Function to generate coordinates (default: nx.spring_layout)
                    Can be nx.spring_layout, nx.circular_layout, nx.random_layout, etc.
    
    Returns:
        Tuple of (modified_graph, list_of_intersections)
    """
    # Auto-detect graph type if not specified
    if is_multigraph is None:
        is_multigraph = isinstance(G, (nx.MultiGraph, nx.MultiDiGraph))
    
    # Check if nodes have coordinates
    has_coords = False
    if G.nodes:
        sample_node = next(iter(G.nodes(data=True)))
        node_data = sample_node[1] if len(sample_node) > 1 else {}
        has_coords = (('x' in node_data and 'y' in node_data) or 
                     ('lon' in node_data and 'lat' in node_data))
        print("coordinates detected")
    
    # Add coordinates if needed
    G_with_coords = G.copy()
    if not has_coords and add_coordinates:
        print(f"Graph has no coordinates. Generating layout...")
        
        # Choose layout function
        if layout_func is None:
            # Default to spring layout for better edge separation
            layout_func = nx.spring_layout
        
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
    
    elif not has_coords and not add_coordinates:
        print("Warning: Graph has no coordinates and add_coordinates=False")
        print("Cannot perform geometric intersection analysis")
        return G.copy(), []
    
    # Convert graph edges to line segments
    segments = []
    
    if is_multigraph:
        # Handle MultiGraph/MultiDiGraph with edge keys
        edge_iter = G_with_coords.edges(keys=True, data=True)
        for edge_data in edge_iter:
            if len(edge_data) == 4:  # (u, v, key, data)
                u, v, key, data = edge_data
                edge_key = (u, v, key)
            else:  # Fallback for unusual cases
                u, v = edge_data[0], edge_data[1]
                key = edge_data[2] if len(edge_data) > 2 else 0
                data = edge_data[-1] if len(edge_data) > 3 else {}
                edge_key = (u, v, key)
            
            if u in G_with_coords.nodes and v in G_with_coords.nodes:
                start_node = G_with_coords.nodes[u]
                end_node = G_with_coords.nodes[v]
                
                # Handle both geographic (lat/lon) and projected coordinates
                if 'y' in start_node and 'x' in start_node:
                    start_point = Point(start_node['x'], start_node['y'])
                    end_point = Point(end_node['x'], end_node['y'])
                elif 'lat' in start_node and 'lon' in start_node:
                    start_point = Point(start_node['lon'], start_node['lat'])
                    end_point = Point(end_node['lon'], end_node['lat'])
                else:
                    # Skip edges without coordinate information
                    continue
                
                segment = Segment(start_point, end_point, edge_key)
                segments.append(segment)
    
    else:
        # Handle regular Graph/DiGraph without edge keys
        for u, v, data in G_with_coords.edges(data=True):
            if u in G_with_coords.nodes and v in G_with_coords.nodes:
                start_node = G_with_coords.nodes[u]
                end_node = G_with_coords.nodes[v]
                
                # Handle both geographic (lat/lon) and projected coordinates
                if 'y' in start_node and 'x' in start_node:
                    start_point = Point(start_node['x'], start_node['y'])
                    end_point = Point(end_node['x'], end_node['y'])
                elif 'lat' in start_node and 'lon' in start_node:
                    start_point = Point(start_node['lon'], start_node['lat'])
                    end_point = Point(end_node['lon'], end_node['lat'])
                else:
                    # Skip edges without coordinate information
                    continue
                
                # For regular graphs, edge key is just (u, v)
                edge_key = (u, v)
                segment = Segment(start_point, end_point, edge_key)
                segments.append(segment)

    nx.draw(G_with_coords, positions, node_size=100)

    print(f"Analyzing {len(segments)} segments for intersections...")
    print(f"Graph type: {'MultiGraph' if is_multigraph else 'Graph'}")
    
    # Find intersections using Bentley-Ottmann algorithm
    intersections = bentley_ottmann_intersections(segments)
    
    print(f"Found {len(intersections)} intersections")
    
    # Create a copy of the graph to modify (preserve coordinates if added)
    G_planar = G_with_coords.copy()
    removed_edges = []
    
    if remove_intersections and intersections:
        # Remove intersecting edges (keep shorter edges preferentially)
        edges_to_remove = set()
    
        for seg1, seg2, intersection_point in intersections:
            # Calculate edge lengths
            len1 = ((seg1.end.x - seg1.start.x)**2 + (seg1.end.y - seg1.start.y)**2)**0.5
            len2 = ((seg2.end.x - seg2.start.x)**2 + (seg2.end.y - seg2.start.y)**2)**0.5
            
            # Remove the longer edge
            edge_to_remove = seg2.edge_key if len1 < len2 else seg1.edge_key
            edges_to_remove.add(edge_to_remove)
        
        # Remove edges from graph based on graph type
        for edge_key in edges_to_remove:
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
        
        print(f"Removed {len(removed_edges)} intersecting edges")
    else:
        
        intersection_list = list(intersections)

        for seg1, seg2, intersection_point in intersection_list:
            node_counter = len(G_planar.nodes)
            G_planar.add_node(f"inter_{node_counter}", x=intersection_point.x, y=intersection_point.y)

            edges_to_remove = set()
            edges_to_remove.add(seg1.edge_key)
            edges_to_remove.add(seg2.edge_key)
        
            # Remove edges from graph based on graph type
            for edge_key in edges_to_remove:
                if is_multigraph:
                    # MultiGraph: edge_key is (u, v, key)
                    u, v, key = edge_key
                    if G_planar.has_edge(u, v, key):
                        G_planar.remove_edge(u, v, key)
                        removed_edges.append(edge_key)
                else:
                    # Regular Graph: edge_key is (u, v)
                    u, v = edge_key
                    print(u, v)
                    if G_planar.has_edge(u, v):
                        G_planar.remove_edge(u, v)
                        removed_edges.append(edge_key)

            points_to_connect = set()
            points_to_connect.add(seg1.edge_key[0])
            points_to_connect.add(seg1.edge_key[1])
            points_to_connect.add(seg2.edge_key[0])
            points_to_connect.add(seg2.edge_key[1])


            for point in points_to_connect:
                G_planar.add_edge(point, f"inter_{node_counter}", key=f"new_{node_counter}")
                if is_multigraph:
                    new_seg = Segment(
                        Point(G_planar.nodes[point]['x'], G_planar.nodes[point]['y']), 
                        Point(G_planar.nodes[f"inter_{node_counter}"]['x'], G_planar.nodes[f"inter_{node_counter}"]['y']),
                        edge_key=(point, f"inter_{node_counter}", f"new_{node_counter}")
                    )
                else:
                    new_seg = Segment(
                        Point(G_planar.nodes[point]['x'], G_planar.nodes[point]['y']), 
                        Point(G_planar.nodes[f"inter_{node_counter}"]['x'], G_planar.nodes[f"inter_{node_counter}"]['y']),
                        edge_key=(point, f"inter_{node_counter}")
                    )
                segments.append(new_seg)
                
                for i in range(0, len(intersection_list)):
                    intersection = intersection_list[i]
                    if intersection == (seg1, seg2, intersection_point):
                        continue
                    must_edit = False
                    for j in range(0, 2):
                        # print(j, intersection[j])
                        if intersection[j].start == new_seg.start or intersection[j].end == new_seg.end or intersection[j].start == new_seg.end or intersection[j].end == new_seg.start:
                            # print("must edit")
                            must_edit = True
                            index = j
                            break
                    if not must_edit: continue

                    # print("editing intersection")
                    if index == 0:
                        intersection_list[i] = (new_seg, intersection[1], intersection[2])
                        # print("New intersection:", intersection_list[i][2])
                    elif index == 1:
                        intersection_list[i] = (intersection[0], new_seg, intersection[2])
                        # print("New intersection:", intersection_list[i][2])

        print(f"Added {node_counter - len(edges_to_remove)} nodes at intersections and connected them")

    # Verify planarity (convert to undirected for planarity check)
    if hasattr(G_planar, 'to_undirected'):
        test_graph = G_planar.to_undirected()
    else:
        test_graph = G_planar
    
    # nx.draw(G_planar, node_size=100)

    is_planar = nx.is_planar(test_graph)
    print(f"Resulting graph is planar: {is_planar}")

    # planar_pos = nx.planar_layout(test_graph)
    # nx.draw(test_graph, planar_pos, node_size=100)
    
    return G_planar, intersections

# Example usage and testing
def test_bentley_ottmann_osmnx():
    """Test the implementation with an OSMnx graph"""
    print("Testing Bentley-Ottmann algorithm with OSMnx graph...")
    
    # Get a small urban area graph
    try:
        point = (40.7831, -73.9712)  # NYC
        G = ox.graph_from_point(point, dist=500, network_type='drive')
        print(f"Original graph: {len(G.nodes)} nodes, {len(G.edges)} edges")
        print(f"Graph type: {type(G).__name__}")
        
        # Test planarity enforcement with auto-detection
        G_planar, intersections = ensure_planar_graph(G, remove_intersections=True)
        print(f"Planar graph: {len(G_planar.nodes)} nodes, {len(G_planar.edges)} edges")
        
        # Verify the result
        is_original_planar = nx.is_planar(G.to_undirected())
        is_result_planar = nx.is_planar(G_planar.to_undirected())
        
        print(f"\nResults:")
        print(f"Original graph planar: {is_original_planar}")
        print(f"Modified graph planar: {is_result_planar}")
        print(f"Intersections found: {len(intersections)}")
        
        return G_planar, intersections
        
    except Exception as e:
        print(f"Error testing with OSMnx: {e}")
        return None, []

def test_different_graph_types():
    """Test with different NetworkX graph types"""
    print("\n" + "="*50)
    print("Testing with different graph types...")
    
    # Create test graphs of different types
    graphs = {
        'Graph': nx.Graph(),
        'DiGraph': nx.DiGraph(), 
        'MultiGraph': nx.MultiGraph(),
        'MultiDiGraph': nx.MultiDiGraph()
    }
    
    # Add some test nodes with coordinates
    test_nodes = [
        (1, {'x': 0, 'y': 0}),
        (2, {'x': 1, 'y': 1}), 
        (3, {'x': 0, 'y': 1}),
        (4, {'x': 1, 'y': 0})
    ]
    
    for graph_name, G in graphs.items():
        print(f"\nTesting {graph_name}:")
        
        # Add nodes
        for node_id, attrs in test_nodes:
            G.add_node(node_id, **attrs)
        
        # Add edges (different syntax for multi vs regular graphs)
        if 'Multi' in graph_name:
            G.add_edge(1, 2, key='edge1')
            G.add_edge(3, 4, key='edge2') 
            G.add_edge(1, 4, key='edge3')
        else:
            G.add_edge(1, 2)
            G.add_edge(3, 4)
            G.add_edge(1, 4)
        
        try:
            # Test with explicit multigraph parameter
            is_multi = 'Multi' in graph_name
            G_planar, intersections = ensure_planar_graph(G, is_multigraph=is_multi)
            
            print(f"  Original edges: {G.number_of_edges()}")
            print(f"  Planar edges: {G_planar.number_of_edges()}")
            print(f"  Intersections: {len(intersections)}")
            
        except Exception as e:
            print(f"  Error: {e}")

if __name__ == "__main__":
    # Run the tests
    planar_graph, intersections = test_bentley_ottmann_osmnx()
    test_different_graph_types()