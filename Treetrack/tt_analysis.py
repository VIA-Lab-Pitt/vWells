from tt_stats import *
from tt_skeleton import *
from VIA_methods import *


class vWellAnalysis:  # For managing the segmentation
    def __init__(self, index_image, vWell_graph, resampling_factor, spacing, input_paths, FromSaveObject):
        # Variables
        self.input_paths = input_paths
        self.resampling_factor = resampling_factor

        # Folders for analysis stuff
        self.input_paths.append(f"{self.input_paths[2]}_ANALYSIS_RF{self.resampling_factor}")

        # More variables
        self.file_extension = None
        self.FromSave = FromSaveObject
        self.resampled_input_image_spacing = spacing
        self.results_handling = vesselResults(FromSaveObject)
        self.vessel_key_counter = 0

        # Analysis / Image stuff
        self.graph_of_vWells = vWell_graph
        self.graph_of_fragments = None
        self.VesselSkeletonizeAlgorithm = None
        self.vWell_index_image = index_image

        self.path_image = None
        self.vWells_in_path = None

    # SKELETONIZE VESSEL ANALYSIS ALGORITHM
    def ManageSkeletonVesselAnalysis(self, show_skeleton, itk_input_image, line_width):
        print(f"Starting vessel analysis...")
        list_of_root_points, skeleton_image = None, None
        start = time.perf_counter_ns()

        # Reset the graph structure and skeletonize object
        self.graph_of_fragments = vesselFragment()
        self.VesselSkeletonizeAlgorithm = vesselSkeletonizeAlgorithm(self.input_paths, self.resampling_factor, self.FromSave, self.graph_of_fragments, self.path_image)

        # Check for saved data
        self.file_extension = "_AIRF"
        list_of_root_points, data_exists = self.VesselSkeletonizeAlgorithm.CheckForSavedData(self.file_extension)

        if not data_exists:
            # Skeletonize the segmentation
            skeleton_image = self.VesselSkeletonizeAlgorithm.SkeletonizeSegment()

            # Find branch and boundary points from skeleton
            list_of_root_points = self.VesselSkeletonizeAlgorithm.FindRootPoints()

        # Get the actors
        if show_skeleton:
            # Marching cubes the skeleton
            skeleton_actor = Get3DVessel(skeleton_image * 2, self.resampled_input_image_spacing, 1)

            return list_of_root_points, [[0, skeleton_actor]], data_exists
        else:
            if not data_exists:
                # Use the direction to find the fragments of vessels
                self.VesselSkeletonizeAlgorithm.FindVesselFragments()

                # Get the not fitted line actors
                branch_actors = GetLineActors(self.graph_of_fragments.fragment_pixel_locations, self.resampling_factor, itk_input_image, line_width)
            else:
                # Get the fitted line actors
                branch_actors = GetLineActors(self.graph_of_fragments.fitted_points, self.resampling_factor, itk_input_image, line_width)

        end = time.perf_counter_ns()
        print(f'\nTotal time taken for vessel analysis was {(end - start) / 1_000_000_000:.3f} seconds.')
        return list_of_root_points, branch_actors, data_exists

    def FinishVesselAnalysis(self, itk_input_image, line_width, r2_threshold):
        print(f"Finishing vessel analysis...")
        start = time.perf_counter_ns()
        # Do the pre-calculations necessary
        self.VesselSkeletonizeAlgorithm.FragmentPreCalculations(r2_threshold, self.resampled_input_image_spacing[0])

        # Save the data
        self.VesselSkeletonizeAlgorithm.SaveData(self.file_extension)

        # Get the branch line actors
        branch_actors = GetLineActors(self.graph_of_fragments.fitted_points, self.resampling_factor, itk_input_image, line_width)

        end = time.perf_counter_ns()
        print(f'\nTotal time taken for finished vessel analysis was {(end - start) / 1_000_000_000:.3f} seconds.')
        return branch_actors

    def SplitFragment(self, new_branch_point, itk_input_image, line_width):
        # Call the specific version
        root_points, self.graph_of_fragments, valid_new_branch_point = self.VesselSkeletonizeAlgorithm.SplitFragmentAndRerunFinder(new_branch_point)

        # Get the skeleton line actors
        branch_actors = GetLineActors(self.graph_of_fragments.fragment_pixel_locations, self.resampling_factor, itk_input_image, line_width)

        return root_points, branch_actors, valid_new_branch_point

    def GetAndSaveFragmentInformation(self, points):
        # Get the individual shortest paths
        mm_per_pixel = self.resampled_input_image_spacing[0]

        total_path_of_fragments = []
        total_path_of_nodes = []
        for i in range(len(points) - 1):
            point1 = points[i]
            point2 = points[i + 1]
            list_of_fragments, shortest_path = self.graph_of_fragments.dijkstra_shortest_path(point1, point2)

            if not list_of_fragments or not shortest_path:
                print(f"No path found between {point1} and {point2}")
                return []

            total_path_of_fragments.extend(list_of_fragments)
            total_path_of_nodes.extend(shortest_path)

        # Check to see if stored info is reversed or not
        fragment_info_reversed = CheckFragmentInfoReversed(total_path_of_fragments, self.graph_of_fragments.fragment_pixel_locations, total_path_of_nodes)

        # Get stats
        total_distance = 0
        total_diameters = []
        total_curvature = []
        for i, fragment in enumerate(total_path_of_fragments):
            total_distance += self.graph_of_fragments.fragment_distances.get(fragment) * mm_per_pixel
            fragment_diameters = self.graph_of_fragments.fragment_diameters.get(fragment).copy()
            fragment_curvature = self.graph_of_fragments.fragment_curvature.get(fragment).copy()
            # Make sure the information is in the order that you picked the points in
            if fragment_info_reversed[i]:
                fragment_diameters.reverse()
                fragment_curvature.reverse()
            # Append to total path list
            if i == 0:
                total_diameters.extend(fragment_diameters)
                total_curvature.extend(fragment_curvature)
            else:
                total_diameters.extend(fragment_diameters[1:])
                total_curvature.extend(fragment_curvature[1:])

        # Summary calculations
        print(f"\nFor the: {self.results_handling.vessels[self.vessel_key_counter]}")
        # If the vessel doesn't exist
        if total_distance == 0:
            total_distance = np.nan
            average_diameter = np.nan
            print(f"\tSkipped...")
        # If the vessel is there
        else:
            average_diameter = np.average(total_diameters) * mm_per_pixel
            print(f"\tmm per pixel = {mm_per_pixel:.3f}")
            print(f"\tAverage diameter of selected fragment is: {average_diameter:.2f} mm")
            print(f"\tLength of selected path is: {total_distance:.2f} mm")
            # PlotEnumeratedValues(total_diameters, "Path Diameter")
            # PlotEnumeratedValues(total_curvature, "Path curvature")

        # Add everything to the data frame
        self.results_handling.AddData(self.vessel_key_counter, "Diameter (mm)", np.array(total_diameters) * mm_per_pixel)
        self.results_handling.AddData(self.vessel_key_counter, "Average Diameter (mm)", average_diameter)
        self.results_handling.AddData(self.vessel_key_counter, "Curvature (1/mm)", total_curvature)
        self.results_handling.AddData(self.vessel_key_counter, "Length", total_distance)
        self.results_handling.AddData(self.vessel_key_counter, "mm per pixel", mm_per_pixel)

        # Next vessel
        self.vessel_key_counter += 1

        return total_path_of_fragments

    def ShowAndSaveVesselData(self):
        # Final calculations comparing the whole brain
        self.results_handling.FinalBrainCalculations()

        # Save everything
        file_extension  = "_RESULTS_RF"
        self.results_handling.SaveResults(self.input_paths, self.resampling_factor, file_extension)

        # Reset in house data frame
        self.results_handling.ResetDataFrame()
        self.vessel_key_counter = 0

    def InteriorVsExteriorVoxelDistributions(self, itk_image):
        print("Plotting interior vs exterior voxel distributions...")
        # Convert ITK image to numpy (resampled if needed)
        if self.resampling_factor == 1:
            intensity_image = itk.array_from_image(itk_image)
        else:
            intensity_image, _ = DoubleSampleITKImage(itk_image)

        # Make sure shapes match
        index_image = self.vWell_index_image
        if index_image.shape != intensity_image.shape:
            raise ValueError("Index image and intensity image must have same dimensions.")

        # Segmentation blob vWells
        blob_vWells = set(self.vWells_in_path)

        # Immediate background vWells
        background_vWells = set()
        for vWell in blob_vWells:
            background_vWells.update(set(self.graph_of_vWells.neighbors_of(vWell)) - blob_vWells)

        # Build mask
        interior_mask = np.isin(index_image, list(blob_vWells))
        exterior_mask = np.isin(index_image, list(background_vWells))

        interior_voxels = intensity_image[interior_mask]
        exterior_voxels = intensity_image[exterior_mask]

        print(f"\tTotal Interior voxels: {interior_voxels.size}, "
              f"Total Exterior voxels: {exterior_voxels.size}")

        return interior_voxels, exterior_voxels

    def GetSegmentedVesselActor(self, opacity):
        # Make an image out of the list of shortest paths
        self.path_image = GetPathImage(self.vWell_index_image, self.vWells_in_path, self.graph_of_vWells.vWell_pixel_locations)

        # Get the marching cubes thing
        segmented_vessel_actor = Get3DVessel(self.path_image, self.resampled_input_image_spacing, opacity)

        return segmented_vessel_actor

    def LoadProgress(self):
        directory_path = self.input_paths[0]
        list_of_paths_extension = "_LOPRF"
        # Load in segmentation
        loaded_arrays, paths_exist = self.FromSave.LoadInSavedNumpyArrays(directory_path, list_of_paths_extension)

        if paths_exist:
            list_of_paths = loaded_arrays['arr_0']
            # Get vWells from list of paths
            list_of_path_vWells = [vWell for path in list_of_paths for vWell in path[1]]
            self.vWells_in_path = set(list_of_path_vWells)
            print(f"\tSuccessfully loaded in saved data.")
            return True
        else:
            return False