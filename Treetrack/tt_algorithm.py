from tt_graphs import vWellGraph
from tt_segment import vWellSegment
from tt_analysis import vWellAnalysis
from VIA_methods import *
import numpy as np
from collections import defaultdict
from numpy.lib.stride_tricks import sliding_window_view
import time
import itk
from numba import njit, prange
from numba.typed import Dict, List
from numba.types import UniTuple, int32, boolean


# Numba optimized functions (must be outside class)
@njit(cache=True)
def fast_compute_direction(variance_image, offsets, pointing_at_you):
    h, w, d = variance_image.shape
    direction_image = np.zeros((h, w, d), dtype=np.int32)

    for i in range(h):
        for j in range(w):
            for k in range(d):
                last_variance = variance_image[i, j, k]
                for direction_index in range(1, len(offsets)):
                    a = i + offsets[direction_index, 0]
                    b = j + offsets[direction_index, 1]
                    c = k + offsets[direction_index, 2]

                    if 0 <= a < h and 0 <= b < w and 0 <= c < d:
                        if direction_image[a, b, c] == pointing_at_you[direction_index]:
                            continue

                        val = variance_image[a, b, c]
                        if val <= last_variance:
                            direction_image[i, j, k] = direction_index
                            last_variance = val

    return direction_image


@njit(cache=True)
def fast_find_vwells(direction_image, minima_coords, offsets, pointing_at_you):
    h, w, d = direction_image.shape
    # Initialize to -1 so that 0 is not confused with empty space
    index_image = np.full((h, w, d), -1, dtype=np.int32)

    # Preallocate stacks for the floodfill
    max_elements = h * w * d
    stack_x = np.zeros(max_elements, dtype=np.int32)
    stack_y = np.zeros(max_elements, dtype=np.int32)
    stack_z = np.zeros(max_elements, dtype=np.int32)

    for minima in range(len(minima_coords)):
        stack_pointer = 0
        stack_x[stack_pointer] = minima_coords[minima, 0]
        stack_y[stack_pointer] = minima_coords[minima, 1]
        stack_z[stack_pointer] = minima_coords[minima, 2]
        stack_pointer += 1

        while stack_pointer > 0:
            stack_pointer -= 1
            a = stack_x[stack_pointer]
            b = stack_y[stack_pointer]
            c = stack_z[stack_pointer]

            # If current pixel is not filled, fill
            if index_image[a, b, c] == -1:
                index_image[a, b, c] = minima
            else:
                continue

            for direction_index in range(1, len(offsets)):
                x = a + offsets[direction_index, 0]
                y = b + offsets[direction_index, 1]
                z = c + offsets[direction_index, 2]

                # Make sure neighbor is within bounds
                if 0 <= x < h and 0 <= y < w and 0 <= z < d:
                    # If neighbor is pointing at you, add to floodfill stack
                    if direction_image[x, y, z] == pointing_at_you[direction_index]:
                        stack_x[stack_pointer] = x
                        stack_y[stack_pointer] = y
                        stack_z[stack_pointer] = z
                        stack_pointer += 1

    return index_image


@njit(cache=True)
def fast_find_neighbors(index_image, kernel_radius, edges_dict):
    height, width, depth = index_image.shape

    # Iterate through every voxel in the 3D image
    for y in range(height):
        for x in range(width):
            for z in range(depth):

                current_id = index_image[y, x, z]

                # Skip empty voxels
                if current_id == -1:
                    continue

                # Calculate the exact bounding box for the neighborhood.
                # Using max() and min() prevents out-of-bounds indexing and
                # eliminates the need for boundary checks inside the loops.
                y_start = max(0, y - kernel_radius)
                y_end = min(height, y + kernel_radius + 1)

                x_start = max(0, x - kernel_radius)
                x_end = min(width, x + kernel_radius + 1)

                z_start = max(0, z - kernel_radius)
                z_end = min(depth, z + kernel_radius + 1)

                # Iterate directly over the valid neighbor coordinates
                for ny in range(y_start, y_end):
                    for nx in range(x_start, x_end):
                        for nz in range(z_start, z_end):

                            neighbor_id = index_image[ny, nx, nz]

                            # Ignore empty voxels and self-references
                            if neighbor_id != -1 and neighbor_id != current_id:

                                # Always format (small, large) to avoid double-counting edges
                                if current_id < neighbor_id:
                                    edges_dict[(current_id, neighbor_id)] = True
                                else:
                                    edges_dict[(neighbor_id, current_id)] = True


# Rebuild vWell_pixel_locations from the saved index image (avoids pickling the dict)
def reconstruct_pixel_locations(vWell_index_image):
    valid_mask = vWell_index_image != -1
    indices = np.argwhere(valid_mask)
    ids = vWell_index_image[valid_mask]

    # Group coordinates by their vWell ID
    sort_order = np.argsort(ids)
    sorted_ids = ids[sort_order]
    sorted_indices = indices[sort_order]
    unique_ids, split_indices = np.unique(sorted_ids, return_index=True)
    grouped_coords = np.split(sorted_indices, split_indices[1:])

    # Store an (N, 3) ndarray per ID - iterable as `for x, y, z in <values>`
    pixel_locations = defaultdict(set)
    for ID, coords in zip(unique_ids, grouped_coords):
        pixel_locations[ID] = coords

    return pixel_locations


class vWellAlgorithm:  # For managing the vWell algorithm
    def __init__(self):
        # Variables
        self.input_paths = None
        self.FromSave = None
        self.numpy_resampled_input_image = None
        self.resampled_input_image_spacing = None
        self.resampling_factor = None
        self.kernel_radius = None
        self.neighbor_offsets = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 1), (0, 0, -1),  # core and centers
                                 (-1, -1, 0), (0, -1, 1), (1, -1, 0), (0, -1, -1),  # U edges
                                 (-1, 0, 1), (1, 0, 1), (1, 0, -1), (-1, 0, -1),  # E edges
                                 (-1, 1, 0), (0, 1, 1), (1, 1, 0), (0, 1, -1),  # D edges
                                 (-1, -1, 1), (1, -1, 1), (1, -1, -1), (-1, -1, -1),  # U corners
                                 (-1, 1, 1), (1, 1, 1), (1, 1, -1), (-1, 1, -1)  # D corners
                                 )
        self.pointing_at_you = (0, 3, 4, 1, 2, 6, 5, 17, 18, 15, 16, 13, 14, 11, 12, 9, 10, 7, 8, 25, 26, 23, 24, 21, 22, 19, 20)

        # Image stuff
        self.Segmentation = None
        self.Analysis = None
        self.graph_of_vWells = vWellGraph()
        self.vWell_index_image = None
        self.list_of_minima = None

    # VWELL ALGORITHM
    def StartVwellFindingAlgorithm(self, input_paths, input_image, resampling_factor, kernel_radius, build_segment_object=False, build_analysis_object=False):
        # Running and managing the vWell algorithm
        start = time.perf_counter_ns()

        # Set the class variables
        self.resampling_factor = resampling_factor
        self.kernel_radius = kernel_radius

        # Data Handling class
        self.input_paths = input_paths
        self.input_paths.append(f"{self.input_paths[2]}_DATA_RF{resampling_factor}")
        self.FromSave = DataHandling(self.input_paths[2], resampling_factor)

        # Set image spacing
        self.resampled_input_image_spacing = np.array(input_image.GetSpacing()) / resampling_factor

        # Check for saved data
        file_extension = "_SERF"
        numpy_variance_image, vWell_mean_image, data_exists = self.CheckForSavedData(file_extension)

        if not data_exists:
            # Resample image according to resampling factor
            if resampling_factor == 2:
                self.numpy_resampled_input_image, self.resampled_input_image_spacing = DoubleSampleITKImage(input_image)
            else:  # Resampling factor = 1
                self.numpy_resampled_input_image = itk.array_from_image(input_image).astype(np.uint32)

            # Load or compute the variance image
            numpy_variance_image = self.ComputeVarianceImage()

            # Find the direction image of the variance image
            direction_image, self.list_of_minima = self.ComputeDirectionImage(numpy_variance_image)

            # Find vWells and fill them
            self.vWell_index_image = self.FindVwells(direction_image)

            # Find vWell neighbors
            self.FindVwellNeighbors()

            # Find vWell mean image and standard deviation image
            vWell_mean_image = self.FindVwellMeanStdevTotal()

            # Save everything necessary for future use
            self.SaveData(file_extension, numpy_variance_image, vWell_mean_image)

        # If we are segmenting make the segmentation object
        if build_segment_object:
            self.BuildSegmentObject()
        elif build_analysis_object:
            self.BuildAnalysisObject()

        end = time.perf_counter_ns()
        print(f'\nTotal time taken for vWell algorithm was {(end - start) / 1_000_000_000:.3f} seconds.\n')

        # Convert back to ITK from numpy
        itk_variance_image = NumpyToITKandSetSpacing(numpy_variance_image, self.resampled_input_image_spacing)
        itk_mean_image = NumpyToITKandSetSpacing(vWell_mean_image, self.resampled_input_image_spacing)

        return itk_variance_image, np.max(numpy_variance_image), itk_mean_image

    def CheckForSavedData(self, file_extension):
        print("Trying to retrieve saved files...")
        # Check to see if data folder exists
        directory_path = os.path.join(self.input_paths[0], self.input_paths[3])
        if not os.path.exists(directory_path):
            print(f"\tSaved data does not exist in: {directory_path}")
            return None, None, False

        # Check to see if loaded data exists
        loaded_images, images_exist = self.FromSave.LoadInSavedNumpyArrays(directory_path, file_extension)

        loaded_dicts, dicts_exist = self.FromSave.LoadInSavedDefaultDicts(directory_path,  file_extension, total_dicts=6)

        if images_exist and dicts_exist:
            numpy_variance_image = loaded_images['arr_0']
            vWell_mean_image = loaded_images['arr_1']
            self.vWell_index_image = loaded_images['arr_2']
            self.graph_of_vWells.graph = loaded_dicts[0]
            self.graph_of_vWells.vWell_means = loaded_dicts[1]
            self.graph_of_vWells.vWell_stdevs = loaded_dicts[2]
            self.graph_of_vWells.vWell_total_pixels = loaded_dicts[3]
            self.graph_of_vWells.vWell_roots = loaded_dicts[4]
            # Slot 5 (vWell_pixel_locations) is reconstructed from the index image
            # instead of loaded from disk - see SaveData
            self._reconstruct_pixel_locations()
            print(f"\tSuccessfully loaded in saved data.")
            return numpy_variance_image, vWell_mean_image, True
        else:
            return None, None, False

    def _reconstruct_pixel_locations(self):
        # Rebuild vWell_pixel_locations from vWell_index_image after a cache load
        print("\tReconstructing vWell_pixel_locations from index image...")
        start = time.perf_counter_ns()
        self.graph_of_vWells.vWell_pixel_locations = reconstruct_pixel_locations(self.vWell_index_image)
        end = time.perf_counter_ns()
        print(f'\tReconstruction took: {(end - start) / 1_000_000_000:.3f} seconds.')

    def SaveData(self, file_extension, numpy_variance_image, vWell_mean_image):
        print("Saving computed data structures...")
        # vWell_pixel_locations is redundant with vWell_index_image and pickling it
        # OOMs on large scans. Save an empty placeholder; reconstruct on load.
        empty_pixel_locations = defaultdict(set)

        graph_structure = [
            self.graph_of_vWells.graph,
            self.graph_of_vWells.vWell_means,
            self.graph_of_vWells.vWell_stdevs,
            self.graph_of_vWells.vWell_total_pixels,
            self.graph_of_vWells.vWell_roots,
            empty_pixel_locations]

        # Combine images and graph arrays into one list
        images = [
            numpy_variance_image,
            vWell_mean_image,
            self.vWell_index_image]

        # Check if the new directory exists
        directory_path = os.path.join(self.input_paths[0], self.input_paths[3])
        if not os.path.exists(directory_path):
            # If not, create the directory
            os.makedirs(directory_path)

        # Save the numpy images together
        self.FromSave.SaveNumpyArrays(directory_path, file_extension, images)

        # Save the default dicts together
        self.FromSave.SaveDefaultDicts(directory_path, file_extension, graph_structure)

        print("\tComputed data saved successfully")

    def ComputeVarianceImage(self):
        print("Computing variance image...")
        start = time.perf_counter_ns()

        # Use float64 to prevent overflow when squaring large voxel values
        image = self.numpy_resampled_input_image.astype(np.float64)
        radius = self.kernel_radius
        kernel_size = (2 * radius) + 1
        neighborhood_size = kernel_size ** 3

        # Pad the image
        # NumPy calls this mode='edge' (repeats the outermost edge values)
        padded_img = np.pad(image, pad_width=radius, mode='edge')

        # Create sliding window views for the image and the squared image
        # Shape becomes (X, Y, Z, k, k, k)
        windows_image = sliding_window_view(padded_img, window_shape=(kernel_size, kernel_size, kernel_size))
        windows_squared_image = sliding_window_view(padded_img**2, window_shape=(kernel_size, kernel_size, kernel_size))

        # Sum over the last 3 axes (which represent our k x k x k neighborhood cube)
        sum_of_image = np.sum(windows_image, axis=(3, 4, 5))
        sum_of_squares = np.sum(windows_squared_image, axis=(3, 4, 5))

        # Calculate integer variance
        # N * sum(x^2) - (sum(x))^2
        variance_image = (neighborhood_size * sum_of_squares) - (sum_of_image**2)

        end = time.perf_counter_ns()
        print(f"\tMax Variance: {np.max(variance_image)}")
        print(f'\tVariance image calculation took: {(end - start) / 1_000_000_000:.3f} seconds.')

        return variance_image

    def ComputeDirectionImage(self, variance_image):
        print("Finding direction image...")
        start = time.perf_counter_ns()

        # Prep variables for Numba
        offsets_array = np.array(self.neighbor_offsets, dtype=np.int32)
        pointing_array = np.array(self.pointing_at_you, dtype=np.int32)

        # Compute direction image using Numba. This was previously slow due to a bunch of nested for loops.
        direction_image = fast_compute_direction(variance_image, offsets_array, pointing_array)

        # Vectorized list of minima
        self.list_of_minima = np.argwhere(direction_image == 0)

        # Count the number of zeros, avoiding redundant checks
        num_zeros = len(self.list_of_minima)
        print(f"\tNumber of zeros in direction image: {num_zeros}")
        print(f"\tLocal minima found = {num_zeros}")

        end = time.perf_counter_ns()
        print(f'\tFinding direction image took: {(end - start) / 1_000_000_000:.3f} seconds.')
        return direction_image, self.list_of_minima

    def FindVwells(self, direction_image):
        print("Finding vWells...")
        start = time.perf_counter_ns()

        offsets_array = np.array(self.neighbor_offsets, dtype=np.int32)
        pointing_array = np.array(self.pointing_at_you, dtype=np.int32)

        # Floodfill with Numba
        vWell_index_image = fast_find_vwells(direction_image, self.list_of_minima, offsets_array, pointing_array)

        # Update root dictionary for graph
        for minima in range(len(self.list_of_minima)):
            self.graph_of_vWells.vWell_roots[minima] = tuple(self.list_of_minima[minima])

        # Build pixel_locations as a dict of ndarrays (much faster than the
        # old per-voxel tuple loop, downstream iteration is unchanged)
        self.graph_of_vWells.vWell_pixel_locations = reconstruct_pixel_locations(vWell_index_image)

        end = time.perf_counter_ns()
        print(f'\tFloodfilling took: {(end - start) / 1_000_000_000:.3f} seconds.')

        return vWell_index_image

    def FindVwellNeighbors(self):
        print("Finding vWell neighbors...")
        start = time.perf_counter_ns()

        # Create the numba dictionary in Python space
        edges_dict = Dict.empty(
            key_type=UniTuple(int32, 2),
            value_type=boolean
        )

        fast_find_neighbors(self.vWell_index_image, self.kernel_radius, edges_dict)

        # Iterate over the dictionary keys
        for node1, node2 in edges_dict.keys():
            self.graph_of_vWells.add_edge(node1, node2)

        end = time.perf_counter_ns()
        print(f'\tFinding vWell neighbors took: {(end - start) / 1_000_000_000:.3f} seconds.')

    def FindVwellMeanStdevTotal(self):
        print("Finding means and standard deviations...")
        start = time.perf_counter_ns()

        # Extract only populated pixels
        ids = self.vWell_index_image.ravel()
        vals = self.numpy_resampled_input_image.ravel()

        valid_mask = ids != -1
        valid_ids = ids[valid_mask]
        valid_vals = vals[valid_mask]

        # Use np.bincount to aggregate data across IDs instantly
        counts = np.bincount(valid_ids)
        sums = np.bincount(valid_ids, weights=valid_vals)
        sums_squares = np.bincount(valid_ids, weights=valid_vals.astype(np.float64) ** 2)

        unique_ids = np.unique(valid_ids)

        # Track means for the mapped 3D output image
        means_map = np.zeros(np.max(unique_ids) + 1, dtype=np.float64)

        # Populate dictionaries in the graph class
        for ID in unique_ids:
            n = counts[ID]
            mean = sums[ID] / n
            variance = (sums_squares[ID] / n) - (mean ** 2)
            stdev = np.sqrt(max(variance, 0))

            self.graph_of_vWells.vWell_means[ID] = mean
            self.graph_of_vWells.vWell_stdevs[ID] = stdev
            self.graph_of_vWells.vWell_total_pixels[ID] = n

            means_map[ID] = mean

        # Rebuild the final 3D mean image via array mapping
        vWell_mean_image = np.where(self.vWell_index_image != -1, means_map[self.vWell_index_image], 0.0)

        end = time.perf_counter_ns()
        print(f'\tFinding means and standard deviations took: {(end - start) / 1_000_000_000:.3f} seconds.')

        return vWell_mean_image

    def BuildSegmentObject(self):
        print("Running tt_segment...")
        self.Segmentation = vWellSegment(self.vWell_index_image, self.graph_of_vWells, self.resampling_factor, self.resampled_input_image_spacing, self.input_paths, self.FromSave)

    def BuildAnalysisObject(self):
        print("Running tt_analysis...")
        self.Analysis = vWellAnalysis(self.vWell_index_image, self.graph_of_vWells, self.resampling_factor, self.resampled_input_image_spacing, self.input_paths, self.FromSave)