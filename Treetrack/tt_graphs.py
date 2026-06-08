# Classes for VIA Lab
from VIA_methods import *
from collections import defaultdict
import numpy as np
from heapq import heappop, heappush
from numba import njit

# These functions are optimized using numba and must be kept outside the class
# Standard T-test, used for shortest path
@njit(cache=True)
def fast_ttest_from_stats(mean1, std1, n1, mean2, std2, n2):
    if n1 <= 1 or n2 <= 1:
        return 0.0

    var1 = std1 ** 2
    var2 = std2 ** 2

    # Denominator for Welch's T-Test, with two weighted variances
    denominator = np.sqrt((var1 / n1) + (var2 / n2))
    if denominator == 0:
        return 0.0

    return np.abs(mean1 - mean2) / denominator


# T-test from running sums, used for region grow and fill from point
@njit(cache=True)
def fast_ttest_from_sums(n1, S1, SS1, n2, S2, SS2):
    if n1 <= 1 or n2 <= 1:
        return 0.0

    mean1 = S1 / n1
    mean2 = S2 / n2

    var1 = (SS1 - (S1 * S1 / n1)) / (n1 - 1)
    var2 = (SS2 - (S2 * S2 / n2)) / (n2 - 1)

    denominator = np.sqrt((var1 / n1) + (var2 / n2))
    if denominator == 0:
        return 0.0

    return np.abs(mean1 - mean2) / denominator


class vWellGraph:  # For graph structures
    def __init__(self):
        self.graph = defaultdict(set)
        self.vWell_means = defaultdict(np.float64)
        self.vWell_stdevs = defaultdict(np.float64)
        self.vWell_total_pixels = defaultdict(int)
        self.vWell_roots = defaultdict(tuple)
        self.vWell_pixel_locations = defaultdict(set)

    def add_edge(self, node1, node2):
        self.graph[node1].add(node2)
        self.graph[node2].add(node1)

    def neighbors_of(self, node):
        return self.graph.get(node).copy()

    # Shortest path stuff
    def calculate_distance(self, node1, node2):
        # distance is the t_value between two nodes
        mean1 = self.vWell_means.get(node1)
        mean2 = self.vWell_means.get(node2)

        stdev1 = self.vWell_stdevs.get(node1)
        stdev2 = self.vWell_stdevs.get(node2)

        n1 = self.vWell_total_pixels.get(node1)
        n2 = self.vWell_total_pixels.get(node2)

        # For some reason, numba doesn't like when you label your parameters (std1=stdev1). Maybe this is because it
        # works using C?
        distance = fast_ttest_from_stats(mean1, stdev1, n1, mean2, stdev2, n2)
        # New fast t-test already does abs value

        return distance

    def distance_heuristic(self, node_a, node_b):
        # Grab the (x, y, z) root coordinates
        pos_a = np.array(self.vWell_roots[node_a])
        pos_b = np.array(self.vWell_roots[node_b])

        # Calculate the straight-line distance between the roots
        return np.linalg.norm(pos_a - pos_b)

    def a_star_shortest_path(self, start_node, end_node, heuristic_weight):
        # Initialize data structures
        open_set = [(0, start_node)]
        current_lowest_t_scores = {node: float('inf') for node in self.graph}
        current_lowest_t_scores[start_node] = 0
        came_from = {}
        visited_nodes = set()

        while open_set:
            current_cost, current_node = heappop(open_set)

            # If we hit the target, stop completely and return the path
            if current_node == end_node:
                break

            if current_node in visited_nodes:
                continue

            visited_nodes.add(current_node)

            for neighbor in self.neighbors_of(current_node):
                if neighbor in visited_nodes:
                    continue

                # t-score to travel to node
                step_cost = self.calculate_distance(current_node, neighbor)
                tentative_t_score = current_lowest_t_scores[current_node] + step_cost

                if tentative_t_score < current_lowest_t_scores.get(neighbor, float('inf')):
                    came_from[neighbor] = current_node
                    current_lowest_t_scores[neighbor] = tentative_t_score
                    heuristic_distance = self.distance_heuristic(neighbor, end_node)

                    # Apply heuristic weight here if needed (may be needed for super blurry scan with super low t_values at boundary, but I have not come across such a scan.)
                    total_cost = tentative_t_score + (heuristic_distance * heuristic_weight)
                    heappush(open_set, (total_cost, neighbor))

        # If the heap empties, we never found a path
        if current_lowest_t_scores[end_node] == float('inf'):
            print(f"\tNo path found")
            return None  # No path found

        # Reconstruct path
        shortest_path = []
        current = end_node
        while current != start_node:
            shortest_path.append(current)
            current = came_from[current]
        shortest_path.append(start_node)
        shortest_path.reverse()

        return shortest_path

    def HemisphereCheck(self, vWell, neighbor_vWells, points, vectors, sphere_radius):
        good_neighbors = set()
        root_in_hemisphere = False
        # Check to see if neighbors are outside of sphere if vWell is in sphere
        for point_number, point in enumerate(points):
            vWell_to_point_distance, vWell_to_point_vector = DistanceAndVectorBetweenPoints(point, self.vWell_roots.get(vWell))

            # Continue to next point if root vWell is outside of sphere (good vWell)
            if vWell_to_point_distance > sphere_radius:
                continue

            # Continue to the next point if the root vWell is NOT hemisphere (good vWell)
            vWell_point_dot = np.dot(vectors[point_number], vWell_to_point_vector)
            if vWell_point_dot > 0:
                continue

            # If the vWell makes it here, it does not pass the hemisphere check
            root_in_hemisphere = True

            # Check the vWell's neighbors since it is in the bad hemisphere
            for neighbor in neighbor_vWells:
                neighbor_to_point_distance, neighbor_to_point_vector = DistanceAndVectorBetweenPoints(point, self.vWell_roots.get(neighbor))
                # # If the neighbor is inside the sphere,
                if neighbor_to_point_distance < sphere_radius:
                    # good_neighbors.add(neighbor)

                    # If the neighbor is in the good hemisphere
                    neighbor_to_point_dot = np.dot(vectors[point_number], neighbor_to_point_vector)
                    if neighbor_to_point_dot > 0:
                        # Then it is a good neighbor
                        good_neighbors.add(neighbor)

        # If the vWell is not in either hemisphere, then all its neighbors are good
        if not root_in_hemisphere:
            good_neighbors.update(neighbor_vWells)

        return good_neighbors

    # Precompute stats for running sums
    def PrecomputeStats(self):
        self.vWell_sums = {}
        self.vWell_sum_sqs_intra = {} # Intra-sample only (Region Grow)
        self.vWell_sum_sqs_full = {}  # Inter and intra-sample variance (Hole Fill)

        for node in self.graph:
            n = self.vWell_total_pixels.get(node, 0)
            mean = self.vWell_means.get(node, 0.0)
            stdev = self.vWell_stdevs.get(node, 0.0)

            # Sum of intensities
            self.vWell_sums[node] = mean * n

            # Sum of squares (intra-sample only)
            self.vWell_sum_sqs_intra[node] = n * (mean ** 2)

            # Sum of squares (inter-sample and intra-sample)
            self.vWell_sum_sqs_full[node] = n * ((stdev ** 2) + (mean ** 2))

    def CompareBlobAndBackground(self, blob_n, blob_S, blob_SS, bg_n, bg_S, bg_SS):
        # Replaced everything with the numba function
        return fast_ttest_from_sums(blob_n, blob_S, blob_SS, bg_n, bg_S, bg_SS)
