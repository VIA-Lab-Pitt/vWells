# Classes for VIA Lab
from VIA_methods import *
from tt_fragments import *
import skimage.morphology
from scipy.ndimage import label
from scipy.ndimage import distance_transform_edt

class vesselSkeletonizeAlgorithm:
    def __init__(self, input_paths, resampling_factor, FromSaveObject, graph_of_fragments, segment):
        # Member variables
        self.input_paths = input_paths
        self.resampling_factor = resampling_factor
        self.FromSave = FromSaveObject
        self.graph_of_fragments = graph_of_fragments
        self.segment = (segment == 2).astype(np.uint8)
        self.root_points = []
        self.skeleton_image = None
        self.skeleton_radii_image = None

        # Direction stuff
        self.neighbor_offsets = ((1, 0, 0), (0, 1, 0), (0, 0, 1),  # FRD centers
                                 (1, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 1),  # FRD corner
                                 (0, 0, -1),  # L center
                                 (0, 1, -1), (1, 0, -1), (1, 1, -1),  # FLD corner
                                 (-1, 0, 0),  # B center
                                 (-1, 1, 0), (-1, 0, 1), (-1, 1, 1),  # RBD corner
                                 (0, -1, 0),  # U center
                                 (1, -1, 0), (0, -1, -1), (1, -1, -1),  # FLU corner
                                 (0, -1, 1), (1, -1, 1),  # FRU corner
                                 (-1, -1, 0), (-1, -1, 1),  # RBU corner
                                 (-1, 0, -1), (-1, 1, -1),  # LBD corner
                                 (-1, -1, -1))  # LBU corner

    def CheckForSavedData(self, file_extension):
        print("Trying to retrieve saved files...")
        # Check to see if data folder exists
        directory_path = os.path.join(self.input_paths[0], self.input_paths[4])
        if not os.path.exists(directory_path):
            print(f"\tSaved data does not exist in: {directory_path}")
            return self.root_points, False

        # Check to see if loaded data exists
        loaded_images, images_exist = self.FromSave.LoadInSavedNumpyArrays(directory_path, file_extension)

        loaded_dicts, dicts_exist = self.FromSave.LoadInSavedDefaultDicts(directory_path, file_extension, total_dicts=8)

        if images_exist and dicts_exist:
            self.skeleton_image = loaded_images['arr_0']
            self.skeleton_radii_image = loaded_images['arr_1']
            root_points = loaded_images['arr_2']
            # Convert to a list of tuples
            self.root_points = [tuple(row) for row in root_points]
            self.graph_of_fragments.roots_neighbors = loaded_dicts[0]
            self.graph_of_fragments.connection_fragments = loaded_dicts[1]
            self.graph_of_fragments.fragment_pixel_locations = loaded_dicts[2]
            self.graph_of_fragments.fragment_end_points = loaded_dicts[3]
            self.graph_of_fragments.fitted_points = loaded_dicts[4]
            self.graph_of_fragments.fragment_distances = loaded_dicts[5]
            self.graph_of_fragments.fragment_diameters = loaded_dicts[6]
            self.graph_of_fragments.fragment_curvature = loaded_dicts[7]

            print(f"\tSuccessfully loaded in old vessel analysis saved data.")
            return self.root_points, True
        else:
            return self.root_points, False

    def SaveData(self, file_extension):
        print("Saving computed data structures...")
        graph_structure = [
            self.graph_of_fragments.roots_neighbors,
            self.graph_of_fragments.connection_fragments,
            self.graph_of_fragments.fragment_pixel_locations,
            self.graph_of_fragments.fragment_end_points,
            self.graph_of_fragments.fitted_points,
            self.graph_of_fragments.fragment_distances,
            self.graph_of_fragments.fragment_diameters,
            self.graph_of_fragments.fragment_curvature]

        # Combine images and graph arrays into one list
        images = [
            self.skeleton_image,
            self.skeleton_radii_image,
            self.root_points]

        # Check if the new directory exists
        directory_path = os.path.join(self.input_paths[0], self.input_paths[4])
        if not os.path.exists(directory_path):
            # If not, create the directory
            os.makedirs(directory_path)

        # Save the numpy images together
        self.FromSave.SaveNumpyArrays(directory_path, file_extension, images)

        # Save the default dicts together
        self.FromSave.SaveDefaultDicts(directory_path, file_extension, graph_structure)

        print("\tComputed old vessel analysis data saved successfully")

    def SkeletonizeSegment(self):
        print("Finding medialness...")
        start = time.perf_counter_ns()

        # Skeletonize path (medialness)
        self.skeleton_image = skimage.morphology.skeletonize(self.segment).astype(np.uint8)

        # Find the skeleton's radii
        # Compute the Euclidean distance transform
        distance_transform = distance_transform_edt(self.segment)

        # Map the distance transform values to the skeleton
        self.skeleton_radii_image = distance_transform * self.skeleton_image

        end = time.perf_counter_ns()
        print(f"\tTotal points in skeleton image: {np.sum(self.skeleton_image == 1)}")
        print(f'\tTotal time taken for skeletonization was {(end - start) / 1_000_000_000:.3f} seconds.')
        return self.skeleton_image

    def FindRootPoints(self):
        print("Finding root points...")
        start = time.perf_counter_ns()
        image_height, image_width, image_depth = self.skeleton_image.shape

        # Kernel size for neighborhood inspection
        kernel_size = 3

        # Pad the image to handle edge cases
        padded_skeleton_image = np.pad(self.skeleton_image, ((1, 1), (1, 1), (1, 1)), mode='constant', constant_values=0)

        # Initialize list to store branch points
        # Go through each voxel in the medialness image
        for i in range(image_height):
            for j in range(image_width):
                for k in range(image_depth):
                    # If the current pixel isn't part of the path
                    if self.skeleton_image[i, j, k] == 0:
                        continue

                    # Extract the local neighborhood
                    neighborhood = np.array(padded_skeleton_image[i:(i + kernel_size), j:(j + kernel_size), k:(k + kernel_size)].copy())

                    # if (i, j, k) == (35, 36, 110):
                    #     print(f"neighborhood: {neighborhood}")
                    #     Visualize3DGrid(neighborhood)

                    # Reset center to 0
                    neighborhood[1, 1, 1] = 0

                    # Use 3D connectivity to find connected components in the neighborhood
                    structure = np.array([[[0, 0, 0],
                                           [0, 1, 0],
                                           [0, 0, 0]],

                                          [[0, 1, 0],
                                           [1, 1, 1],
                                           [0, 1, 0]],

                                          [[0, 0, 0],
                                           [0, 1, 0],
                                           [0, 0, 0]]])

                    # Connected components
                    _, num_features = label(neighborhood, structure=structure)
                    # Any branch points in neighborhood
                    # neighboring_boundary_points = np.sum(neighborhood == 2)
                    neighboring_branch_points = np.sum(neighborhood == 3)

                    # This means not a branch point
                    if num_features == 2:
                        continue

                    # This is a boundary point, allowed always
                    if num_features == 1:
                        padded_skeleton_image[i + 1, j + 1, k + 1] = 2
                        self.root_points.append((i, j, k))

                    # This is a branch point, not allowed next to another branch point
                    # if num_features > 2 and neighboring_branch_points <= 1:
                    #     padded_skeleton_image[i + 1, j + 1, k + 1] = 3
                    #     root_points.append((i, j, k))

                    # This is a branch point
                    if num_features > 2:
                        padded_skeleton_image[i + 1, j + 1, k + 1] = 3
                        self.root_points.append((i, j, k))

        end = time.perf_counter_ns()
        print(f"\tTotal root points found: {len(self.root_points)}")
        print(f'\tTotal time taken for finding branch points was {(end - start) / 1_000_000_000:.3f} seconds.')
        return self.root_points

    def FindVesselFragments(self):
        print("Finding vessel fragments...")
        start = time.perf_counter_ns()
        image_height, image_width, image_depth = self.skeleton_image.shape

        set_of_roots = set(self.root_points)
        unvisited_roots = set(self.root_points)
        root_stack = [self.root_points[0]]
        visited_pixels = set()

        # Process fragments starting from each root in the stack.
        while root_stack:
            start_pixel = root_stack.pop()
            if start_pixel in unvisited_roots:
                unvisited_roots.remove(start_pixel)
            current_pixel = start_pixel
            vector_of_pixels_in_fragment = []
            fragment_complete = False
            end_pixel = None
            adjacent_roots = set()
            counter = 0

            # Iterating through the skeleton, building each fragment
            while True:
                counter += 1
                a, b, c = current_pixel
                found_next_pixel = False
                next_pixel = None

                # Look at each neighbor of the current pixel.
                for offset in self.neighbor_offsets:
                    x = a + offset[0]
                    y = b + offset[1]
                    z = c + offset[2]
                    neighbor_coordinate = (x, y, z)

                    # Skip if we've already visited this pixel.
                    if neighbor_coordinate in visited_pixels:
                        continue

                    # Check if neighbor is within the image bounds.
                    if not (0 <= x < image_height and 0 <= y < image_width and 0 <= z < image_depth):
                        continue

                    # If the neighbor is a root
                    if neighbor_coordinate in set_of_roots:
                        # If you are at the start root, and you find adjacent roots and skip them
                        if counter == 1:
                            adjacent_roots.add(neighbor_coordinate)
                            continue
                        # If you are not at the start root, and you find a root that isn't one of the original adjacent roots
                        elif neighbor_coordinate not in adjacent_roots:
                            # Pick it immediately...
                            next_pixel = neighbor_coordinate
                            found_next_pixel = True
                            break  # Prioritize the root neighbor.
                        # If you make it here, you are probably at the second step where you have found an adjacent root to bounce too, so we skip
                        else:
                            continue

                    # Otherwise, if it is part of the skeleton, record it if none was chosen yet.
                    if self.skeleton_image[x, y, z] == 1 and not found_next_pixel:
                        next_pixel = neighbor_coordinate
                        found_next_pixel = True
                        # Do not break; we want to keep looking in case a root is found.

                # Record the current pixel as part of the fragment.
                vector_of_pixels_in_fragment.append(current_pixel)
                visited_pixels.add(current_pixel)

                # If no neighbor was found, consider the fragment complete using the current pixel.
                if not found_next_pixel:
                    if counter == 1 and adjacent_roots:
                        end_pixel = adjacent_roots.pop()
                        # We are done floodfilling and don't want to look at any more branches coming off of it for now
                        vector_of_pixels_in_fragment.append(end_pixel)
                        fragment_complete = True
                    else:
                        end_pixel = current_pixel
                        # fragment_complete = True
                        fragment_complete = False
                    break

                # If the chosen neighbor is a root, we treat that as the fragment’s end.
                if next_pixel in set_of_roots:
                    end_pixel = next_pixel

                    # Count how many branches come off from this root.
                    remaining_branches = 0
                    e, f, g = end_pixel
                    for offset in self.neighbor_offsets:
                        x = e + offset[0]
                        y = f + offset[1]
                        z = g + offset[2]

                        # Check if we have visited the pixel yet
                        if (x, y, z) in visited_pixels:
                            continue
                        if 0 <= x < image_height and 0 <= y < image_width and 0 <= z < image_depth:
                            # Count the total number of features
                            if self.skeleton_image[x, y, z] == 1:
                                remaining_branches += 1

                    if remaining_branches:
                        # Add as many copies of this root as there are branches.
                        root_stack.extend([end_pixel] * remaining_branches)
                    elif end_pixel in root_stack:
                        # If this root already exists in the stack from another branch, remove the duplicate.
                        root_stack.remove(end_pixel)
                    # We are done floodfilling
                    vector_of_pixels_in_fragment.append(end_pixel)
                    fragment_complete = True
                    break
                else:
                    # Continue along the fragment.
                    current_pixel = next_pixel

            # If the fragment is complete and has more than one pixel, add it to the graph.
            if fragment_complete:
                fragment_ID = self.graph_of_fragments.add_connection_between(start_pixel, end_pixel)
                self.graph_of_fragments.fragment_pixel_locations[fragment_ID].extend(vector_of_pixels_in_fragment)
                self.graph_of_fragments.fragment_end_points[fragment_ID].update({start_pixel, end_pixel})

            # Remove all roots from the visited set to avoid blocking valid paths.
            visited_pixels.difference_update(set_of_roots)

            # If there are still unvisited roots and nothing left in the stack, add another.
            if unvisited_roots and not root_stack:
                root_stack.append(unvisited_roots.pop())

        end = time.perf_counter_ns()
        print(f"\tTotal points in fragmentation: {len(visited_pixels) + len(set_of_roots)}")
        print(f"\tTotal vessel fragments found: {len(self.graph_of_fragments.fragment_pixel_locations)}")
        print(f'\tFinding vessel fragments took: {(end - start) / 1_000_000_000:.3f} seconds.')

    def CombineConnectedFragments(self, combine_fragments, itk_input_image, line_width):
        # Get path of nodes through the fragments
        start_nodes = self.graph_of_fragments.fragment_end_points.get(combine_fragments[0]).copy()
        total_path_of_nodes = list(start_nodes)
        for i in range(len(combine_fragments) - 1):
            next_fragment = combine_fragments[i + 1]
            first_node = set(total_path_of_nodes[-2:])
            next_nodes = self.graph_of_fragments.fragment_end_points.get(next_fragment).copy()
            arranged_list = ArrangeTwoSets3Elements(first_node, next_nodes)
            if not arranged_list:
                print(f"\tNo connection found")
                return False
            else:
                total_path_of_nodes[-2:] = arranged_list

        # Check to see if stored info is reversed or not
        fragment_info_reversed = CheckFragmentInfoReversed(combine_fragments, self.graph_of_fragments.fragment_pixel_locations, total_path_of_nodes)

        # Get combined data
        total_pixel_locations = []
        for i, fragment in enumerate(combine_fragments):
            fragment_pixels = self.graph_of_fragments.fragment_pixel_locations.get(fragment).copy()
            # Make sure the information is in the order that you picked the points in
            if fragment_info_reversed[i]:
                fragment_pixels.reverse()
            # Append to total path list
            if i == 0:
                total_pixel_locations.extend(fragment_pixels)
            else:
                total_pixel_locations.extend(fragment_pixels[1:])

        # Remove old connections
        for i in range(len(total_path_of_nodes) - 1):
            node1 = total_path_of_nodes[i]
            node2 = total_path_of_nodes[i + 1]
            delete_fragment = self.graph_of_fragments.get_connection_between(node1, node2)
            self.graph_of_fragments.remove_connection_between(node1, node2)
            del self.graph_of_fragments.fragment_pixel_locations[delete_fragment]
            del self.graph_of_fragments.fragment_end_points[delete_fragment]

        # Add new connection
        node1 = total_path_of_nodes[0]
        node2 = total_path_of_nodes[-1]
        fragment_ID = self.graph_of_fragments.add_connection_between(node1, node2)
        self.graph_of_fragments.fragment_pixel_locations[fragment_ID].extend(total_pixel_locations)
        self.graph_of_fragments.fragment_end_points[fragment_ID].update({node1, node2})

        # Get the skeleton line actors
        branch_actors = GetLineActors(self.graph_of_fragments.fragment_pixel_locations, self.resampling_factor, itk_input_image, line_width)

        return branch_actors

    def SplitFragmentAndRerunFinder(self, new_branch_point):
        # Check to see if proposed branch point is on skeleton
        x, y, z = new_branch_point
        if not self.skeleton_image[x, y, z]:
            return self.root_points, self.graph_of_fragments, False

        # If it is on the skeleton, add it as a new branch point
        self.root_points.append(new_branch_point)

        # Then rerun the fragment finding
        self.graph_of_fragments = vesselFragment()  # Reset the graph structure
        self.FindVesselFragments()

        return self.root_points, self.graph_of_fragments, True

    def FragmentPreCalculations(self, r2_threshold, mm_per_pixel):
        print("Doing fragment pre-calculations...")
        start = time.perf_counter_ns()

        # Fill in the fragment information dictionaries
        for fragment, pixels_in_fragment in self.graph_of_fragments.fragment_pixel_locations.items():
            list_of_diameters = []
            # Get the fitted line points
            total_fitted_points = 100
            max_polynomial_degree = 9
            fitted_points, coefficients, t_bounds = Fit3DPolynomial(pixels_in_fragment, max_polynomial_degree, total_fitted_points, r2_threshold)
            total_arc_length = ArcLengthFromFittedPolynomial(coefficients, t_bounds)
            curvature_values = CurvatureFromFittedPolynomial(coefficients, t_bounds, 40, mm_per_pixel)

            # Go through each pixel
            for pixel in range(len(pixels_in_fragment)):
                i, j, k = pixels_in_fragment[pixel]
                # Get diameter
                list_of_diameters.append(self.skeleton_radii_image[i, j, k] * 2)

            # Update dictionaries structure
            self.graph_of_fragments.polynomial_coefficients[fragment].extend(coefficients)
            self.graph_of_fragments.fragment_distances[fragment] = total_arc_length
            self.graph_of_fragments.fitted_points[fragment].extend(fitted_points)
            self.graph_of_fragments.fragment_diameters[fragment].extend(list_of_diameters)
            self.graph_of_fragments.fragment_curvature[fragment].extend(curvature_values)

        end = time.perf_counter_ns()
        print(f'\tFragment pre-calculations took: {(end - start) / 1_000_000_000:.3f} seconds.')

    def VisualizeSkeletonAndRadii(self):
        # Get the coordinates and y values of the skeleton points
        skeleton_coords = np.argwhere(self.skeleton_image)
        skeleton_radii = self.skeleton_radii_image[self.skeleton_image > 0]

        Show3DPlot(skeleton_coords, skeleton_radii)
