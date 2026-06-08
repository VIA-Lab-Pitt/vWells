from tt_algorithm import *
from VIA_methods import *

class TreetrackAnalysisGUI:
    def __init__(self):
        self.current_image = None
        self.current_plane_orientation = None
        self.plane_widget_on = True
        self.have_vessel_actor = False
        self.displaying_vessel_actor = False
        self.displaying_sphere_actors = False
        self.displaying_line_actors = False
        self.completed_vWell_algorithm = False
        self.completed_skeleton_analysis = False
        self.completed_vessel_calculations = False
        self.reset_save_fragment = False
        self.reset_combine_fragments = False

        self.sphere_actors = []
        self.sphere_actor_points = []
        self.picked_nodes = []
        self.combine_fragments = []
        self.line_actors = None
        self.path_fragments = []
        self.segmented_vessel_actor = None
        self.vessel_opacity = 0.5
        self.total_screenshots = 0

    def KeypressCommands(self, key_sym):
        # Settings
        resampling_factor = 2
        vWell_kernel_radius = 1
        line_width = 7
        r2_threshold = 0.99
        show_skeleton = False   # True for marching cubes skeleton, false for smoothed skeleton
        # Must use smoothed skeleton for further analysis...

        if key_sym == "Up":
            # Plane widget slice up
            new_slice_index = plane_widget.GetSliceIndex() + 1
            plane_widget.SetSliceIndex(new_slice_index)
            render_window.Render()

        elif key_sym == "Down":
            # Plane widget slice down
            new_slice_index = plane_widget.GetSliceIndex() - 1
            plane_widget.SetSliceIndex(new_slice_index)
            render_window.Render()

        elif key_sym == "Right":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Go to the next image
            self.current_image = IncrementImageListIndex(plane_widget_images, self.current_image, up=True)
            # Update the window and level accordingly
            UpdatePlaneWidgetImageWindowLevel(plane_widget, plane_widget_images, self.current_image)
            # Set the new slice index from stored values * resampling factor
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])

            render_window.Render()

        elif key_sym == "Left":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Go to the previous image
            self.current_image = IncrementImageListIndex(plane_widget_images, self.current_image, down=True)
            # Update the window and level accordingly
            UpdatePlaneWidgetImageWindowLevel(plane_widget, plane_widget_images, self.current_image)
            # Set the new slice index from stored values * resampling factor
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])

            render_window.Render()

        elif key_sym == "a" and self.line_actors:
            # Toggle on and off the line actors
            self.displaying_line_actors = not self.displaying_line_actors
            for line in self.line_actors:
                line[1].SetVisibility(self.displaying_line_actors)
                line[1].Modified()
            render_window.Render()

        elif key_sym == "A" and self.sphere_actors:
            # Toggle on and off the sphere actor
            self.displaying_sphere_actors = not self.displaying_sphere_actors
            for sphere in self.sphere_actors:
                sphere.SetVisibility(self.displaying_sphere_actors)
                sphere.Modified()
            render_window.Render()

        elif key_sym == "b" and self.completed_skeleton_analysis and not show_skeleton:
            if self.reset_save_fragment:
                # Reset lines
                for fragment in self.path_fragments:
                    index = next((index for index, (ID, actor) in enumerate(self.line_actors) if ID == fragment), None)
                    line_actor = self.line_actors[index][1]
                    line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(GetVTKColor(fragment)))
                self.path_fragments = []
                # Reset points
                for point in range(len(self.picked_nodes)):
                    # Essentially the undo for b
                    RemoveLastPointResetActor(self.picked_nodes, self.sphere_actor_points, self.sphere_actors, "Red")
                self.reset_save_fragment = False

            # Pick the actor
            picked_actor = PickCurrentActor(render_window_interactor, renderer)
            # For spheres
            if picked_actor in self.sphere_actors and self.completed_vessel_calculations:
                picked_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Blue"))
                actor_index = self.sphere_actors.index(picked_actor)
                # Append the numpy root point corresponding to the sphere
                numpy_root_point = self.sphere_actor_points[actor_index]
                self.picked_nodes.append(numpy_root_point)

            # For lines
            if not self.completed_vessel_calculations:
                actor_fragment_ID = next((fragment_ID for fragment_ID, actor in self.line_actors if actor == picked_actor), None)
                print(f"actor_fragment_ID: {actor_fragment_ID}")
                if actor_fragment_ID is not None:
                    picked_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Blue"))
                    self.combine_fragments.append(actor_fragment_ID)

        elif key_sym == "B" and self.completed_skeleton_analysis:
            # Essentially the undo for b
            if self.picked_nodes and self.completed_vessel_calculations:
                RemoveLastPointResetActor(self.picked_nodes, self.sphere_actor_points, self.sphere_actors, "Red")

            if self.combine_fragments and not self.completed_vessel_calculations:
                RemoveLastLineResetActor(self.combine_fragments, self.line_actors)

        elif key_sym == "C" and self.completed_skeleton_analysis and not show_skeleton:
            # Only do anything if a line actor is selected yk
            if self.combine_fragments:
                # Remove the old the lines
                for actor in self.line_actors:
                    renderer.RemoveActor(actor[1])

                # Combine the selected fragments and update necessary things
                line_actors = vWell_algorithm.Analysis.VesselSkeletonizeAlgorithm.CombineConnectedFragments(self.combine_fragments, itk_input_image, line_width)
                if not line_actors:
                    for fragment in self.combine_fragments:
                        index = next((index for index, (ID, actor) in enumerate(self.line_actors) if ID == fragment), None)
                        line_actor = self.line_actors[index][1]
                        line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(GetVTKColor(fragment)))
                else:
                    self.line_actors = line_actors

                # Reset list of fragments to combine
                self.combine_fragments = []

                # Display the lines
                for actor in self.line_actors:
                    renderer.AddActor(actor[1])
                self.displaying_line_actors = True

        elif key_sym == "d" and self.completed_vWell_algorithm and self.have_vessel_actor:
            # Plot segmentation vs background histogram
            interior_voxels, exterior_voxels = vWell_algorithm.Analysis.InteriorVsExteriorVoxelDistributions(itk_input_image)

            # Plot the histogram
            bins = 200
            plt.figure(figsize=(8, 5))
            plt.hist(
                interior_voxels,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=2,
                label="Interior voxels"
            )
            plt.hist(
                exterior_voxels,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=2,
                label="Exterior voxels"
            )

            plt.xlabel("Voxel intensity")
            plt.ylabel("Probability density")
            plt.title("Interior vs Exterior Voxel Intensity Distributions")
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            # plt.savefig(os.path.join(input_paths[0], "blob_background_histogram.png"), dpi=300, bbox_inches="tight")

            print(
                f"\tInterior mean={np.mean(interior_voxels):.3f}, "
                f"std={np.std(interior_voxels):.3f}")
            print(
                f"\tExterior mean={np.mean(exterior_voxels):.3f}, "
                f"std={np.std(exterior_voxels):.3f}")

            plt.show()

        elif key_sym == "h" and self.completed_skeleton_analysis and not show_skeleton:
            # Pick the actor
            picked_actor = PickCurrentActor(render_window_interactor, renderer)
            # For spheres
            if picked_actor in self.sphere_actors:
                actor_index = self.sphere_actors.index(picked_actor)
                # Remove the actor
                renderer.RemoveActor(self.sphere_actors[actor_index])

        elif key_sym == "l" and self.completed_vWell_algorithm:
            print("Trying to load saved progress...")
            # Load in list of paths to continue off from a previous point
            if vWell_algorithm.Analysis.LoadProgress():
                # Replace the vessel actor
                self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Analysis, self.vessel_opacity)
                self.have_vessel_actor = True
                self.displaying_vessel_actor = True

        elif key_sym == "m" and self.completed_vWell_algorithm and self.have_vessel_actor and not self.completed_vessel_calculations:
            # Remove old spheres and lines
            if self.displaying_sphere_actors:
                for actor in self.sphere_actors:
                    renderer.RemoveActor(actor)
                self.sphere_actors = []

            if self.displaying_line_actors:
                for actor in self.line_actors:
                    renderer.RemoveActor(actor[1])
                self.line_actors = []

            # Start the analysis algorithm
            root_points, self.line_actors, data_exists = vWell_algorithm.Analysis.ManageSkeletonVesselAnalysis(show_skeleton, itk_input_image, line_width)
            self.completed_skeleton_analysis = True
            if data_exists:
                self.completed_vessel_calculations = True

            # Add a sphere actor for each branch point
            for point in root_points:
                vtk_point = ConvertNumpyPointToVTK(point, itk_input_image, resampling_factor)
                # Create a sphere at the selected point
                CreateSphereActor(vtk_point, "Red", 0.3, self.sphere_actors, renderer)
                self.sphere_actor_points.append(point)

            self.displaying_sphere_actors = True

            # Display the lines
            for actor in self.line_actors:
                renderer.AddActor(actor[1])
            self.displaying_line_actors = True

            # Toggle on and off the plane widget
            if self.plane_widget_on:
                plane_widget.Off()
            self.plane_widget_on = False

        elif key_sym == "n" and self.completed_skeleton_analysis and not self.completed_vessel_calculations and not show_skeleton:
            # Remove the old lines
            for actor in self.line_actors:
                renderer.RemoveActor(actor[1])

            # Finish the analysis algorithm and get the branch actors
            self.line_actors = vWell_algorithm.Analysis.FinishVesselAnalysis(itk_input_image, line_width, r2_threshold)
            self.completed_vessel_calculations = True

            # Display the fitted lines
            for actor in self.line_actors:
                renderer.AddActor(actor[1])
            self.displaying_line_actors = True

        elif key_sym == "o" and self.have_vessel_actor:
            # Switch between 0.5 opacity and 1
            if self.vessel_opacity == 0.5:
                self.vessel_opacity = 1
            elif self.vessel_opacity == 1:
                self.vessel_opacity = 0.5

            print(f"Vessel actor opacity set to {self.vessel_opacity}")

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Analysis, self.vessel_opacity)
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

        elif key_sym == "P":
            # File path
            ss_file_path = os.path.expanduser(os.path.join(input_paths[0], f"{input_paths[2]}_tt_Analyze_screenshot_{self.total_screenshots}.png"))

            # Temporarily multiply line actor width by 2
            for actor in self.line_actors:
                actor[1].GetProperty().SetLineWidth(line_width * 2)

            render_window.Render()

            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(render_window)
            windowToImageFilter.SetScale(2)
            windowToImageFilter.SetInputBufferTypeToRGB()  # key change
            windowToImageFilter.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(ss_file_path)
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()

            # Reset line widths
            for actor in self.line_actors:
                actor[1].GetProperty().SetLineWidth(line_width)

            self.total_screenshots += 1
            print(f"Screenshot saved to {ss_file_path}")

        elif key_sym == "D":
            # Toggle on and off the plane widget
            if self.plane_widget_on:
                plane_widget.Off()
            elif not self.plane_widget_on:
                plane_widget.On()

            self.plane_widget_on = not self.plane_widget_on

        elif key_sym == "s" and self.completed_vessel_calculations and not show_skeleton:
            # Save the selected path information
            self.path_fragments = vWell_algorithm.Analysis.GetAndSaveFragmentInformation(self.picked_nodes)

            for fragment in self.path_fragments:
                index = next((index for index, (ID, actor) in enumerate(self.line_actors) if ID == fragment), None)
                line_actor = self.line_actors[index][1]
                line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Blue"))
            self.reset_save_fragment = True

        elif key_sym == "S" and self.completed_vessel_calculations and not show_skeleton:
            # Save big brain dictionary and Excel file and reset data frame
            vWell_algorithm.Analysis.ShowAndSaveVesselData()

        elif key_sym == "t" and self.completed_skeleton_analysis and not self.completed_vessel_calculations and not show_skeleton:
            # Use VTKPointPicker to find the nearest skeleton point
            itk_point, vtk_point, numpy_point = GetNearestActorToMouseAsITKVTKNumpyCoordinates(render_window_interactor, renderer, itk_input_image, resampling_factor)

            # Split the fragment by adding a new branch point and re-fragment finding
            root_points, temporary_line_actors, valid_new_branch_point = vWell_algorithm.Analysis.SplitFragment(numpy_point, itk_input_image, line_width)

            if valid_new_branch_point:
                # If it is a valid new branch point, remove old spheres and lines
                if self.displaying_sphere_actors:
                    for actor in self.sphere_actors:
                        renderer.RemoveActor(actor)
                    self.sphere_actors = []

                if self.displaying_line_actors:
                    for actor in self.line_actors:
                        renderer.RemoveActor(actor[1])
                    self.line_actors = []
                # Reset the line actors

                self.line_actors = temporary_line_actors

                # Add a sphere actor for each branch point
                for point in root_points:
                    vtk_point = ConvertNumpyPointToVTK(point, itk_input_image, resampling_factor)
                    # Create a sphere at the selected point
                    CreateSphereActor(vtk_point, "Red", 0.3, self.sphere_actors, renderer)
                    self.sphere_actor_points.append(point)

                self.displaying_sphere_actors = True

                # Display the lines
                for actor in self.line_actors:
                    renderer.AddActor(actor[1])
                self.displaying_line_actors = True

                # Toggle on and off the plane widget
                if self.plane_widget_on:
                    plane_widget.Off()
                self.plane_widget_on = False

        elif key_sym == "u":
            ExportSceneToGLTF(render_window, self.segmented_vessel_actor, self.line_actors, self.sphere_actors, plane_widget)

        elif key_sym == "v" and not self.completed_vWell_algorithm:
            # Start the vWell algorithm (this is the long-running task)
            itk_variance_image, max_variance, itk_mean_image = vWell_algorithm.StartVwellFindingAlgorithm(input_paths, itk_input_image, resampling_factor, vWell_kernel_radius, build_analysis_object=True)

            # Mean image
            flipped_vtk_mean_image = ITKtoVTKandFlipY(itk_mean_image)
            # Append mean image to list of images and data
            plane_widget_list_mean_image = [flipped_vtk_mean_image, 230, 95, resampling_factor]
            plane_widget_images.append(plane_widget_list_mean_image)

            # Variance image
            flipped_vtk_variance_image = ITKtoVTKandFlipY(itk_variance_image)
            # Append variance image to list of images and data
            plane_widget_list_variance_image = [flipped_vtk_variance_image, 0.09 * max_variance, 0.05 * max_variance, resampling_factor]
            plane_widget_images.append(plane_widget_list_variance_image)

            # Update boolean
            self.completed_vWell_algorithm = True

        elif key_sym == "V" and self.have_vessel_actor:
            # Toggle on and off the vessel actor
            self.displaying_vessel_actor = not self.displaying_vessel_actor
            self.segmented_vessel_actor.SetVisibility(self.displaying_vessel_actor)
            self.segmented_vessel_actor.Modified()

        elif key_sym == "x":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to x-axis
            plane_widget.SetPlaneOrientationToXAxes()
            self.current_plane_orientation = 0
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])

        elif key_sym == "y":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to y-axis
            plane_widget.SetPlaneOrientationToYAxes()
            self.current_plane_orientation = 1
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])

        elif key_sym == "z":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to x-axis
            plane_widget.SetPlaneOrientationToZAxes()
            self.current_plane_orientation = 2
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])

        elif key_sym == "question":
            print("Up Arrow    - Increase slice index")
            print("Down Arrow  - Decrease slice index")
            print("Right Arrow - Go to next image")
            print("Left Arrow  - Go to previous image")
            print("a           - Toggle line actors on/off")
            print("A           - Toggle sphere actors on/off")
            print("b           - Select nodes/fragments")
            print("B           - Deselect last node/fragments")
            print("C           - Combine selected fragments")
            print("D           - Toggle plane widget on/off")
            print("h           - Delete sphere actor")
            print("l           - Load in saved segment")
            print("m           - Start analysis algorithm")
            print("n           - Fit polynomials and save fragmentation")
            print("o           - Toggle vessel actor opacity")
            print("P           - Screenshot screen")
            print("s           - Save stats for selected path")
            print("S           - Save all data for vessels")
            print("t           - Split fragment at selected point")
            print("u           - Save current renderer window as OBJ")
            print("v           - Start vWell algorithm")
            print("V           - Toggle 3D vessel actor on/off")
            print("x           - Display slice along x axis")
            print("y           - Display slice along y axis")
            print("z           - Display slice along z axis")
            print("?           - Print this menu")
            print("\n")

        # Always render
        render_window.Render()

    def KeypressCallback(self, caller, ev):
        iren = caller
        key_sym = iren.GetKeySym()
        self.KeypressCommands(key_sym)


# Initialize everything
if __name__ == "__main__":
    # Get the folder and such
    input_paths = GetFileAndFolderNames()
    # input_paths = ['C:/Users/satya/Documents/PycharmProjects/VIA Lab/Test Brain', 'ExtractedCOW.nii', 'ExtractedCOW']
    # Read NIfTI file using ITK
    filename = os.path.join(input_paths[0], input_paths[1])
    itk_input_image = itk.imread(filename)
    print(f"Working on brain scan: {input_paths[2]}")

    # Print size
    input_image_size = itk_input_image.GetLargestPossibleRegion().GetSize()
    print(f"\tImage Dimensions: {input_image_size[0]} x {input_image_size[1]} x {input_image_size[2]}")

    # Reset origin to (0,0,0) and direction matrix to identity
    ResetOriginAndDirection(itk_input_image)

    # Create the vWell segment algorithm
    vWell_algorithm = vWellAlgorithm()

    # Convert ITK image to VTK and flip along y-axis
    flipped_vtk_input_image = ITKtoVTKandFlipY(itk_input_image)

    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(575, 510)

    # Create renderer
    renderer = vtk.vtkRenderer()
    render_window.AddRenderer(renderer)

    # Create interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Create interactor style
    trackball_style = vtk.vtkInteractorStyleTrackballCamera()
    trackball_style.SetCurrentRenderer(renderer)
    render_window_interactor.SetInteractorStyle(trackball_style)

    # Create a list of images and window/level/resampling factor information
    plane_widget_list_original_image = [flipped_vtk_input_image, 230, 95, 1]
    plane_widget_images = [plane_widget_list_original_image]

    # Initialize ttaGUI class
    ttaGUI = TreetrackAnalysisGUI()

    # Set initial plane orientation and image
    ttaGUI.current_plane_orientation = 2
    ttaGUI.current_image = 0

    # Create plane widget
    plane_widget = vtk.vtkImagePlaneWidget()
    plane_widget.SetInteractor(render_window_interactor)
    plane_widget.SetTextureInterpolate(False)
    plane_widget.SetResliceInterpolateToNearestNeighbour()
    plane_widget.SetInputData(plane_widget_images[ttaGUI.current_image][0])  # Initial Image
    plane_widget.SetPlaneOrientationToZAxes()
    plane_widget.SetMarginSizeX(0.0)
    plane_widget.SetMarginSizeY(0.0)

    # Set initial window and level
    plane_widget.SetWindowLevel(plane_widget_images[ttaGUI.current_image][1], plane_widget_images[ttaGUI.current_image][2])

    # Set initial plane index
    vtk_dimensions = flipped_vtk_input_image.GetDimensions()
    image_plane_indices = [vtk_dimensions[0] // 2, vtk_dimensions[1] // 2, vtk_dimensions[2] // 2]
    plane_widget.SetSliceIndex(image_plane_indices[ttaGUI.current_plane_orientation])

    # Initialize the plane widget
    plane_widget.On()

    # Set up callbacks
    render_window_interactor.AddObserver(vtk.vtkCommand.KeyPressEvent, ttaGUI.KeypressCallback)

    # Start the rendering
    render_window.Render()
    render_window_interactor.Initialize()
    render_window_interactor.Start()