# Classes for VIA Lab
from collections import defaultdict
import numpy as np
from heapq import heappop, heappush

class vesselFragment:  # For graph structures for vessels
    def __init__(self):
        # For combining fragments
        self.max_fragments = -1
        # Graph of branch points as nodes
        self.roots_neighbors = defaultdict(set)
        self.connection_fragments = defaultdict(int)
        # Fragment information, keyed by fragment ID
        self.fragment_pixel_locations = defaultdict(list)
        self.fragment_end_points = defaultdict(set)
        self.fitted_points = defaultdict(list)
        self.polynomial_coefficients = defaultdict(list)
        # Stats, keyed by fragment ID
        self.fragment_distances = defaultdict(np.float64)
        self.fragment_diameters = defaultdict(list)
        self.fragment_curvature = defaultdict(list)

    def add_connection_between(self, node1, node2):
        self.max_fragments += 1
        # Use frozenset for the key
        self.connection_fragments[frozenset([node1, node2])] = self.max_fragments
        self.roots_neighbors[node1].add(node2)
        self.roots_neighbors[node2].add(node1)

        return self.max_fragments

    def get_connection_between(self, node1, node2):
        # Use frozenset for the key
        return self.connection_fragments.get(frozenset([node1, node2]))

    def connections_of(self, node):
        return self.roots_neighbors.get(node)

    def remove_connection_between(self, node1, node2):
        # Remove connections
        del self.connection_fragments[frozenset([node1, node2])]
        self.roots_neighbors[node1].remove(node2)
        self.roots_neighbors[node2].remove(node1)

    def dijkstra_shortest_path(self, start_node, end_node):
        # Initialize data structures
        min_heap = [(0, start_node)]
        shortest_distances = {node: float('inf') for node in self.roots_neighbors}
        shortest_distances[start_node] = 0
        last_nodes = defaultdict(int)
        visited_nodes = set()

        while min_heap:
            current_distance, current_node = heappop(min_heap)

            if current_node == end_node:
                break

            if current_distance > shortest_distances.get(current_node):
                continue

            for neighbor in self.connections_of(current_node):
                if neighbor in visited_nodes:
                    continue

                # Get actual distance
                fragment_connection = self.get_connection_between(current_node, neighbor)
                new_distance = self.fragment_distances.get(fragment_connection)

                total_distance = current_distance + new_distance

                if total_distance < shortest_distances.get(neighbor):
                    shortest_distances[neighbor] = total_distance
                    last_nodes[neighbor] = current_node
                    heappush(min_heap, (total_distance, neighbor))

            visited_nodes.add(current_node)

        if shortest_distances[end_node] == float('inf'):
            return [], []  # No path found

        # Reconstruct path
        shortest_path = []
        current = end_node
        while current != start_node:
            shortest_path.append(current)
            current = last_nodes[current]
        shortest_path.append(start_node)
        shortest_path.reverse()

        # Get the fragments from the shortest path
        path_fragments = []
        for i in range(len(shortest_path) - 1):
            node1 = shortest_path[i]
            node2 = shortest_path[i + 1]
            path_fragments.append(self.get_connection_between(node1, node2))

        return path_fragments, shortest_path
