from VIA_methods import *
import numpy as np
from heapq import heappop, heappush
import time


class vWellSegment:  # For managing the segmentation
    def __init__(self, index_image, vWell_graph, resampling_factor, spacing, input_paths, FromSaveObject):
        # Variables
        self.input_paths = input_paths
        self.FromSave = FromSaveObject
        self.resampled_input_image_spacing = spacing
        self.resampling_factor = resampling_factor
        self.sphere_radius = None

        # Image stuff
        self.graph_of_vWells = vWell_graph
        self.vWell_index_image = index_image
        self.list_of_paths = []
        
        # Precompute running sums for region grow
        self.graph_of_vWells.PrecomputeStats()

    # DIJKSTRAS AND REGION GROW
    def DijkstraAndRegionGrow(self, points, directional_bias, max_negative_iterations, percent_drop, sphere_radius):
        # Update variables
        self.sphere_radius = sphere_radius * self.resampling_factor

        # Find the shortest path between two selected paths
        shortest_path = self.FindShortestPath(points, directional_bias)

        # Region grow from the shortest path
        vWells_in_vessel, vWells_around_vessel, t_values, region_filled_vWells = self.RegionGrowFromPath(shortest_path, points, max_negative_iterations, percent_drop)

        # Append the shortest path and final blob to list of paths
        # self.list_of_paths.append([0, vWells_in_vessel])  # 0 for segment
        region_grow_peak_index = len(vWells_in_vessel) - len(shortest_path)
        self.list_of_paths.append([0, set(shortest_path)])  # 0 for shortest path (point based)
        self.list_of_paths.append([2, set(region_filled_vWells[0:region_grow_peak_index])])  # 2 for region grow post-processing

        return t_values, region_filled_vWells, region_grow_peak_index

    def FindShortestPath(self, points, directional_bias):
        print("Finding shortest path...")
        start = time.perf_counter_ns()
        # Find the starting and ending vWells (nodes)
        start_vWell = GetVwellFromPoint(points[0], self.vWell_index_image)
        end_vWell = GetVwellFromPoint(points[1], self.vWell_index_image)

        print(f"\tStart: {start_vWell}, End: {end_vWell}")
        individual_shortest_path = self.graph_of_vWells.a_star_shortest_path(start_vWell, end_vWell, directional_bias)

        end = time.perf_counter_ns()
        print(f"\tvWells in shortest path: {len(individual_shortest_path)}")
        print(f'\tFinding shortest path took: {(end - start) / 1_000_000_000:.3f} seconds.')
        return individual_shortest_path

    def GetBestBackgroundVwell(self, blob, background, points, vectors, blob_stats, background_stats, ss_dict):
        # Running totals
        blob_n, blob_S, blob_SS = blob_stats
        background_n, background_S, background_SS = background_stats
        # Sort the vector by which vWell is best compared to the current background and blob
        stack_heap = []  # min heap to hold indices of local minima
        candidates = list(background)

        for vWell in candidates:
            # Get stats for the candidate
            candidate_n = self.graph_of_vWells.vWell_total_pixels[vWell]
            candidate_S = self.graph_of_vWells.vWell_sums[vWell]
            candidate_SS = ss_dict[vWell]

            # Find new neighbors to add to background
            raw_neighbors = self.graph_of_vWells.neighbors_of(vWell)
            if len(points) > 1:
                valid_neighbors = self.graph_of_vWells.HemisphereCheck(
                    vWell, raw_neighbors.copy(), points, vectors, self.sphere_radius
                )
            else:
                valid_neighbors = raw_neighbors

            neighbors_to_add = valid_neighbors.difference(blob)
            neighbors_to_add.difference_update(background)

            # Get total stats for the new neighbors being added to background
            neighbor_n = sum(self.graph_of_vWells.vWell_total_pixels[n] for n in neighbors_to_add)
            neighbor_S = sum(self.graph_of_vWells.vWell_sums[n] for n in neighbors_to_add)
            neighbor_SS = sum(ss_dict[n] for n in neighbors_to_add)

            # Calculate temporary running totals for the T-Test
            temp_blob_n = blob_n + candidate_n
            temp_blob_S = blob_S + candidate_S
            temp_blob_SS = blob_SS + candidate_SS

            temp_background_n = background_n - candidate_n + neighbor_n
            temp_background_S = background_S - candidate_S + neighbor_S
            temp_background_SS = background_SS - candidate_SS + neighbor_SS

            # T-test
            t_value = self.graph_of_vWells.CompareBlobAndBackground(temp_blob_n, temp_blob_S, temp_blob_SS, temp_background_n,
                                                                    temp_background_S, temp_background_SS)
            # Add the negative t-value and vWell to the heap
            # so that when we heappop, we get the vWell that gives
            # us the largest t-value.
            heappush(stack_heap, (-1 * t_value, vWell))

        if not stack_heap:
            return None
        
        t_value, best_vWell = heappop(stack_heap)

        return best_vWell

    def RegionGrowFromPath(self, shortest_path, points, max_negative_iterations, percent_drop, use_only_intrasample_variance=True):
        print("Region growing...")
        start = time.perf_counter_ns()

        # Decide which variance to use depending on if we are growing from path or point.
        if use_only_intrasample_variance:
            ss_dict = self.graph_of_vWells.vWell_sum_sqs_intra
        else:
            ss_dict = self.graph_of_vWells.vWell_sum_sqs_full

        if len(points) > 1:
            _, vector1 = DistanceAndVectorBetweenPoints(points[0], points[1])
            _, vector2 = DistanceAndVectorBetweenPoints(points[1], points[0])
            vectors = [vector1, vector2]
        else:
            vectors = False

        # Make a set for all the vWells in the shortest path (blob)
        blob = set(shortest_path)

        # Make a set of all the neighbors of the blob (potential background)
        background = set()
        for vWell in blob:
            background.update(self.graph_of_vWells.neighbors_of(vWell))
        # Remove all the vWells that are in the blob
        background.difference_update(blob)
        print(f"\tOriginal vWells in blob: {len(blob)}")
        print(f"\tOriginal vWells in background: {len(background)}")

        # Make a copy in case no neighbors are better:
        return_blob, return_background = blob.copy(), background.copy()

        # Start the region growing
        negative_iterations = 0
        greatest_t_value = 0
        t_values = []
        added_vWells = []

        # blob stats
        blob_n = sum(self.graph_of_vWells.vWell_total_pixels[v] for v in blob)
        blob_S = sum(self.graph_of_vWells.vWell_sums[v] for v in blob)
        blob_SS = sum(ss_dict[v] for v in blob)
        # background stats
        background_n = sum(self.graph_of_vWells.vWell_total_pixels[v] for v in background)
        background_S = sum(self.graph_of_vWells.vWell_sums[v] for v in background)
        background_SS = sum(ss_dict[v] for v in background)


        while True:
            # Compare current blob and background
            anchor_t_value = self.graph_of_vWells.CompareBlobAndBackground(blob_n, blob_S, blob_SS, background_n, background_S, background_SS)

            # Get the best vWell from the current background
            blob_stats = (blob_n, blob_S, blob_SS)
            background_stats = (background_n, background_S, background_SS)
            vWell = self.GetBestBackgroundVwell(blob, background, points, vectors, blob_stats, background_stats, ss_dict)
            added_vWells.append(vWell)

            # Grab stats for the chosen vWell before transferring
            chosen_n = self.graph_of_vWells.vWell_total_pixels[vWell]
            chosen_S = self.graph_of_vWells.vWell_sums[vWell]
            chosen_SS = ss_dict[vWell]

            # Transfer the vWell from background to blob
            if vWell not in blob:
                blob.add(vWell)
            if vWell in background:
                background.remove(vWell)

            # Get neighbors of vWell
            neighbor_vWells = self.graph_of_vWells.neighbors_of(vWell)
            # Remove the neighbors that are already in the blob
            neighbor_vWells.difference_update(blob)
            if len(points) > 1:
                good_neighbors = self.graph_of_vWells.HemisphereCheck(vWell, neighbor_vWells, points, vectors, self.sphere_radius)
                # Find neighboring vWells
                new_neighbors_for_background = good_neighbors.difference(background)
            else:
                new_neighbors_for_background = neighbor_vWells.difference(background)

            # Add only new neighbors to the background
            background.update(new_neighbors_for_background)

            # Get total stats for just the new neighbors being added to background
            neighbor_n = sum(self.graph_of_vWells.vWell_total_pixels[n] for n in new_neighbors_for_background)
            neighbor_S = sum(self.graph_of_vWells.vWell_sums[n] for n in new_neighbors_for_background)
            neighbor_SS = sum(ss_dict[n] for n in new_neighbors_for_background) # Updated dict call

            # Update running totals
            blob_n += chosen_n
            blob_S += chosen_S
            blob_SS += chosen_SS

            background_n = background_n - chosen_n + neighbor_n
            background_S = background_S - chosen_S + neighbor_S
            background_SS = background_SS - chosen_SS + neighbor_SS

            # Compare new background and blob
            new_t_value = self.graph_of_vWells.CompareBlobAndBackground(blob_n, blob_S, blob_SS, background_n, background_S, background_SS)

            # Handle negative iterations
            if new_t_value > anchor_t_value:
                negative_iterations = 0
                # Choosing between the largest maximum or last maximum
                if new_t_value > greatest_t_value:
                    return_blob, return_background = blob.copy(), background.copy()
                    greatest_t_value = new_t_value
            elif new_t_value <= anchor_t_value:
                negative_iterations += 1
            # Implementing both
            if negative_iterations > max_negative_iterations:
                break
            elif new_t_value < percent_drop * greatest_t_value:
                break
            else:
                t_values.append(new_t_value)

        end = time.perf_counter_ns()
        print(f'\tNew length of blob: {len(return_blob)}')
        print(f'\tNew length of background: {len(return_background)}')
        print(f'\tTime taken for region growing was {(end - start) / 1_000_000_000:.3f} seconds.')

        return return_blob, return_background, t_values, added_vWells


    def RegionGrowPointFill(self, point, sphere_radius, max_negative_iterations=0, percent_drop=1):
        # Update variables
        self.sphere_radius = sphere_radius * self.resampling_factor

        # Get starting vWell
        start_vWell = GetVwellFromPoint(point[0], self.vWell_index_image)

        # Region grow from the shortest path, with both inter-sample and intra-sample variance
        vWells_in_vessel, vWells_around_vessel, t_values, region_filled_vWells = self.RegionGrowFromPath([start_vWell], point, max_negative_iterations, percent_drop, use_only_intrasample_variance=False)

        # Append total blob to list of paths
        region_grow_peak_index = len(vWells_in_vessel) - 1
        # self.list_of_paths.append([1, vWells_in_vessel])  # 1 for hole
        self.list_of_paths.append([3, set(region_filled_vWells[0:region_grow_peak_index])])  # 3 for region grow point fill post-processing

        return region_filled_vWells, region_grow_peak_index

    def FillVwell(self, point):
        # Get vWell
        vWell = GetVwellFromPoint(point, self.vWell_index_image)

        # Append vWell to list of paths
        self.list_of_paths.append([1, {vWell}])  # 1 for singular vWell hole

    def FillInnerHoles(self):
        # Get vWells in path from list of paths
        vWells_in_path = set(vWell for path in self.list_of_paths for vWell in path[1])

        # Iterate through all vWells and fill inner holes
        filled_inner_vWells = set()
        corresponding_roots_points = []
        for vWell in vWells_in_path:
            for neighbor in self.graph_of_vWells.neighbors_of(vWell):
                # Check if all second neighbors of the neighbor are in vWells_in_path
                if all(second_neighbor in vWells_in_path for second_neighbor in
                       self.graph_of_vWells.neighbors_of(neighbor)):
                    if neighbor not in vWells_in_path:
                        filled_inner_vWells.add(neighbor)
                        corresponding_roots_points.append(self.graph_of_vWells.vWell_roots.get(neighbor))

        self.list_of_paths.append([1, filled_inner_vWells])  # 1 for holes
        print(f"Filled {len(filled_inner_vWells)} inner holes.")
        return corresponding_roots_points

    def FillAllExceptFrontier(self, start_point):
        # Flood-fill from known vWell in the image frontier. Fill everything else as vessel.
        start_vWell = GetVwellFromPoint(start_point, self.vWell_index_image)
        filled_vWells = set(vWell for path in self.list_of_paths for vWell in path[1])

        # Check if user-selected vWell is already in blob
        if start_vWell in filled_vWells:
            print("Selected point lies inside the filled blob. Please choose a background point.")
            return True

        # Flood-fill all reachable non-blob vWells from start vWell
        frontier = set()
        visited = set()
        to_visit = [start_vWell]

        while to_visit:
            current = to_visit.pop()
            if current in visited or current in filled_vWells:
                continue
            visited.add(current)
            frontier.add(current)

            for neighbor in self.graph_of_vWells.neighbors_of(current):
                if neighbor not in visited and neighbor not in filled_vWells:
                    to_visit.append(neighbor)

        # Segment everything else
        all_vWells = set(self.graph_of_vWells.graph.keys())
        to_fill = all_vWells - frontier
        new_filled_vWells = to_fill - filled_vWells

        # Append vWell to list of paths
        self.list_of_paths.append([1, new_filled_vWells])  # 1 for filled vWell hole

        print(f"Flood-filled background from vWell {start_vWell}.")
        print(f"Filled {len(new_filled_vWells)} vWells.")
        return False

    def RemoveVwell(self, point):
        # Get vWell
        vWell = GetVwellFromPoint(point, self.vWell_index_image)

        # Remove vWells to list of paths
        for path in self.list_of_paths:
            if vWell in path[1]:
                path[1].remove(vWell)

    def RemoveFloatingVwells(self):
        # Get vWells in path from list of paths
        vWells_in_path = set(vWell for path in self.list_of_paths for vWell in path[1])

        removed_vWell_counter = 0
        for vWell in vWells_in_path:
            # If any of the neighbors are in the path: skip
            if any(neighbor in vWells_in_path for neighbor in
                   self.graph_of_vWells.neighbors_of(vWell)):
                continue
            # If not, remove the vWell
            else:
                removed_vWell_counter += 1
                for path in self.list_of_paths:
                    if vWell in path[1]:
                        path[1].remove(vWell)

        print(f"Removed {removed_vWell_counter} floating vWells.")

    def GetSegmentedVesselActor(self, opacity):
        # Get vWells in path from list of paths
        vWells_in_path = set([vWell for path in self.list_of_paths for vWell in path[1]])

        # Make an image out of the list of shortest paths
        path_image = GetPathImage(self.vWell_index_image, vWells_in_path, self.graph_of_vWells.vWell_pixel_locations)

        # Get the marching cubes thing
        segmented_vessel_actor = Get3DVessel(path_image, self.resampled_input_image_spacing, opacity)

        return segmented_vessel_actor

    def SaveLoadProgress(self, save=False, load=False):
        directory_path = self.input_paths[0]
        list_of_paths_extension = "_LOPRF"
        if save:
            # Save the numpy images together
            self.FromSave.SaveNumpyArrays(directory_path, list_of_paths_extension, [np.array(self.list_of_paths)])
        if load:
            loaded_arrays, paths_exist = self.FromSave.LoadInSavedNumpyArrays(directory_path, list_of_paths_extension)

            if paths_exist:
                self.list_of_paths = list(loaded_arrays['arr_0'])
                return True
            else:
                return False
