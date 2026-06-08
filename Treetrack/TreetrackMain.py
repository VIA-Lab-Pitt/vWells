from tt_algorithm import *
from VIA_methods import *

class TreetrackGUI:
    def __init__(self):
        self.current_image = None
        self.current_plane_orientation = None
        self.plane_widget_on = True
        self.have_vessel_actor = False
        self.displaying_vessel_actor = False
        self.displaying_sphere_actor = False
        self.completed_vWell_algorithm = False
        self.just_did_n_or_h = False

        self.region_filled_vWells = []
        self.added_filled_vWells = None
        self.type_of_last_path = None
        self.t_values = None
        self.points = []
        self.sphere_actors = []  # List to store sphere actors
        self.segmented_vessel_actor = None
        self.vessel_opacity = 0.5

    def KeypressCommands(self, key_sym):
        # Settings
        max_negative_iterations = 8
        percent_drop = 0.80
        resampling_factor = 2
        kernel_radius = 1
        aim_intensity = 1
        sphere_radius = 9

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

        elif key_sym == "A" and self.sphere_actors:
            # Toggle on and off the vessel actor
            self.displaying_sphere_actor = not self.displaying_sphere_actor
            for sphere in self.sphere_actors:
                sphere.SetVisibility(self.displaying_sphere_actor)
                sphere.Modified()
            render_window.Render()

        elif key_sym == "b" and self.completed_vWell_algorithm:
            '''Floodfill the entire 'frontier,' which includes everything in the image except for the segmented blob and
            its interior holes, and then add the remaining vWells to the blob.'''

            # Get 'frontier point' from user mouse cursor
            itk_point, vtk_point, numpy_point = GetMousePositionAsITKVTKNumpyCoordinates(
                render_window_interactor, renderer, itk_input_image, resampling_factor)

            # Run the frontier-based filling
            vessel_point = vWell_algorithm.Segmentation.FillAllExceptFrontier(numpy_point)

            if not vessel_point:
                # Update vessel actor
                self.segmented_vessel_actor = ReplaceVesselActor(
                    renderer, self.segmented_vessel_actor,
                    self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)

                # Visualize the chosen 'frontier' point
                CreateSphereActor(vtk_point, "Carrot", 0.22, self.sphere_actors, renderer)

                # Set the flags
                self.displaying_sphere_actor = True
                self.just_did_n_or_h = False
                self.have_vessel_actor = True
                self.displaying_vessel_actor = True

                # Always render
                render_window.Render()

        elif key_sym == "d" and self.completed_vWell_algorithm:
            # Remove the clicked vWell

            # Get vWell under the mouse cursor
            itk_point, vtk_point, numpy_point = GetMousePositionAsITKVTKNumpyCoordinates(render_window_interactor, renderer, itk_input_image, resampling_factor)
            seed_vWell = GetVwellFromPoint(numpy_point, vWell_algorithm.Segmentation.vWell_index_image)

            # Collect the set of currently-segmented vWells (union across all paths)
            segmented_vWells = set()
            for path in vWell_algorithm.Segmentation.list_of_paths:
                segmented_vWells.update(path[1])

            # If the clicked vWell isn't in the segmentation, do nothing
            if seed_vWell not in segmented_vWells:
                print(f"Clicked vWell {seed_vWell} is not in the segmentation — nothing to remove.")
            else:
                # Remove the clicked vWell from every path that contains it
                for path in vWell_algorithm.Segmentation.list_of_paths:
                    path[1].discard(seed_vWell)

                print(f"Removed vWell {seed_vWell}.")

                # Replace the vessel actor
                self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
                self.have_vessel_actor = True
                self.displaying_vessel_actor = True

                # Always render
                render_window.Render()

        elif key_sym == "D":
            # Toggle on and off the plane widget
            if self.plane_widget_on:
                plane_widget.Off()
            elif not self.plane_widget_on:
                plane_widget.On()

            self.plane_widget_on = not self.plane_widget_on

        elif key_sym == "g" and self.have_vessel_actor and self.just_did_n_or_h:
            # Local variable
            list_of_paths = vWell_algorithm.Segmentation.list_of_paths
            if self.added_filled_vWells < len(self.region_filled_vWells):
                if list_of_paths[-1][0] == 2 or list_of_paths[-1][0] == 3:
                    # Remove last path
                    self.type_of_last_path = list_of_paths[-1][0]
                    list_of_paths.pop()

                # Add the next vWell to the path and rerender
                self.added_filled_vWells += 1
                print(f"Showing region filled vWells: {self.added_filled_vWells} of {len(self.region_filled_vWells)}")
                list_of_paths.append([self.type_of_last_path, set(self.region_filled_vWells[0:self.added_filled_vWells])])

            else:
                print(f"Reached maximum region filled vWells")

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "G" and self.have_vessel_actor and self.just_did_n_or_h:
            # Local variable
            list_of_paths = vWell_algorithm.Segmentation.list_of_paths
            if self.added_filled_vWells > 0:
                if list_of_paths[-1][0] == 2 or list_of_paths[-1][0] == 3:
                    # Remove last path
                    self.type_of_last_path = list_of_paths[-1][0]
                    list_of_paths.pop()

                # Add the next vWell to the path and rerender
                self.added_filled_vWells -= 1
                print(f"Showing region filled vWells: {self.added_filled_vWells} of {len(self.region_filled_vWells)}")
                list_of_paths.append([self.type_of_last_path, set(self.region_filled_vWells[0:self.added_filled_vWells])])
            else:
                print(f"Reached minimum region filled vWells")

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "h" and self.completed_vWell_algorithm:
            print(f"Region filling hole from singular point...")
            # Do a mini region fill from a singular point
            itk_point, vtk_point, numpy_point = GetMousePositionAsITKVTKNumpyCoordinates(render_window_interactor, renderer, itk_input_image, resampling_factor)

            # Create a sphere at the selected point
            CreateSphereActor(vtk_point, "Blue", 0.22, self.sphere_actors, renderer)
            self.displaying_sphere_actor = True

            # Region grow from that point
            self.region_filled_vWells, self.added_filled_vWells = vWell_algorithm.Segmentation.RegionGrowPointFill([numpy_point], sphere_radius)

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)

            # Set the flags
            self.just_did_n_or_h = True
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "i" and self.completed_vWell_algorithm:
            # Do a mini region fill from a singular point
            itk_point, vtk_point, numpy_point = GetMousePositionAsITKVTKNumpyCoordinates(render_window_interactor, renderer, itk_input_image, resampling_factor)
            # Create a sphere at the selected point
            CreateSphereActor(vtk_point, "Green", 0.22, self.sphere_actors, renderer)
            self.displaying_sphere_actor = True

            # Fill that singular point
            vWell_algorithm.Segmentation.FillVwell(numpy_point)

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.just_did_n_or_h = False
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            print(f"Filled singular vWell at numpy point: {numpy_point}")

            # Always render
            render_window.Render()

        elif key_sym == "I" and self.completed_vWell_algorithm:
            # Fill inner holes if they exist
            numpy_points = vWell_algorithm.Segmentation.FillInnerHoles()

            # Delete floating vWells if they exist
            vWell_algorithm.Segmentation.RemoveFloatingVwells()

            for point in numpy_points:
                vtk_point = ConvertNumpyPointToVTK(point, itk_input_image, resampling_factor)
                # Create a sphere for the selected points
                CreateSphereActor(vtk_point, "Green", 0.22, self.sphere_actors, renderer)
                self.displaying_sphere_actor = True

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.just_did_n_or_h = False
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "l" and self.completed_vWell_algorithm:
            print("Trying to load saved progress...")
            # Load in list of paths to continue off from a previous point
            if vWell_algorithm.Segmentation.SaveLoadProgress(load=True):
                # Replace the vessel actor
                self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
                self.just_did_n_or_h = False
                self.have_vessel_actor = True
                self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "n" and self.completed_vWell_algorithm:
            itk_point, vtk_point, numpy_point = GetMousePositionAsITKVTKNumpyCoordinates(render_window_interactor, renderer, itk_input_image, resampling_factor)
            print(f"\tITK point selected: {[round(coord, 2) for coord in itk_point]}")
            self.points.append(numpy_point)

            # Create a sphere at the selected point
            CreateSphereActor(vtk_point, "Red", 0.22, self.sphere_actors, renderer)
            self.displaying_sphere_actor = True

            if len(self.points) >= 2:
                # Get the last two points for dijkstras
                dijkstra_points = self.points[-2:]

                # Do the dijkstras
                self.t_values, self.region_filled_vWells, self.added_filled_vWells = vWell_algorithm.Segmentation.DijkstraAndRegionGrow(dijkstra_points, aim_intensity, max_negative_iterations, percent_drop, sphere_radius)

                # Replace the vessel actor
                self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)

                # Set the flags
                self.just_did_n_or_h = True
                self.have_vessel_actor = True
                self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

            # if len(self.points) >= 2:
            #     PlotEnumeratedValues(self.t_values, "t values against indices")

        elif key_sym == "N":
            # Reset the point picking to start a new vessel
            print("Reset points for segmentation.")
            self.just_did_n_or_h = False
            self.points = []

        elif key_sym == "o" and self.have_vessel_actor:
            # Switch between 0.5 opacity and 1
            if self.vessel_opacity == 0.5:
                self.vessel_opacity = 1
            elif self.vessel_opacity == 1:
                self.vessel_opacity = 0.5

            print(f"Vessel actor opacity set to {self.vessel_opacity}")

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # Always render
            render_window.Render()

        elif key_sym == "S" and self.have_vessel_actor:
            # Save list of paths to continue off from a previous point
            vWell_algorithm.Segmentation.SaveLoadProgress(save=True)
            print("Progress Saved.")

        elif key_sym == "Z" and self.have_vessel_actor:
            # Undo the last segment
            print("Removing previous path...")
            if vWell_algorithm.Segmentation.list_of_paths:
                # Remove last path
                removed_path = vWell_algorithm.Segmentation.list_of_paths.pop()
                self.type_of_last_path = removed_path[0]
                # If the path came from dijkstras, remove point
                if removed_path[0] == 0 and self.points:
                    self.points.pop()
                # Check to see if we should turn 'g' off
                if removed_path[0] != 2:
                    self.just_did_n_or_h = False
                # Remove the last sphere
                if self.sphere_actors:
                    if removed_path[0] in (0, 1, 3):
                        last_sphere = self.sphere_actors.pop()
                        renderer.RemoveActor(last_sphere)

            # Replace the vessel actor
            self.segmented_vessel_actor = ReplaceVesselActor(renderer, self.segmented_vessel_actor, self.have_vessel_actor, vWell_algorithm.Segmentation, self.vessel_opacity)
            self.added_filled_vWells = 0
            self.have_vessel_actor = True
            self.displaying_vessel_actor = True

            # If the list of paths is empty
            if not vWell_algorithm.Segmentation.list_of_paths:
                self.have_vessel_actor = False

            # Always render
            render_window.Render()

        elif key_sym == "v" and not self.completed_vWell_algorithm:
            # Start the vWell algorithm (this is the long-running task)
            itk_variance_image, max_variance, itk_mean_image = vWell_algorithm.StartVwellFindingAlgorithm(input_paths, itk_input_image, resampling_factor, kernel_radius, build_segment_object=True)

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

            render_window.Render()

            # Update boolean
            self.completed_vWell_algorithm = True

        elif key_sym == "V" and self.have_vessel_actor:
            # Toggle on and off the vessel actor
            self.displaying_vessel_actor = not self.displaying_vessel_actor
            self.segmented_vessel_actor.SetVisibility(self.displaying_vessel_actor)
            self.segmented_vessel_actor.Modified()
            render_window.Render()

        elif key_sym == "x":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to x-axis
            plane_widget.SetPlaneOrientationToXAxes()
            self.current_plane_orientation = 0
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])
            render_window.Render()

        elif key_sym == "y":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to y-axis
            plane_widget.SetPlaneOrientationToYAxes()
            self.current_plane_orientation = 1
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])
            render_window.Render()

        elif key_sym == "z":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to x-axis
            plane_widget.SetPlaneOrientationToZAxes()
            self.current_plane_orientation = 2
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])
            render_window.Render()

        elif key_sym == "question":
            print("Up Arrow    - Increase slice index")
            print("Down Arrow  - Decrease slice index")
            print("Right Arrow - Go to next image")
            print("Left Arrow  - Go to previous image")
            print("A           - Toggle sphere actors on/off")
            print("d           - Remove clicked vWell")
            print("D           - Toggle plane widget on/off")
            print("g           - Add singular vWells from region grow")
            print("G           - Remove singular vWells from region grow")
            print("h           - Mini region fill for holes")
            print("i           - Fill in singular vWell")
            print("I           - Fill inner holes")
            print("l           - Load in saved segment")
            print("n           - Pick continuous points to connect")
            print("N           - Start a new path to connect")
            print("o           - Toggle vessel actor opacity")
            print("S           - Save entire current segment")
            print("v           - Start vWell algorithm")
            print("V           - Toggle 3D vessel actor on/off")
            print("x           - Display slice along x axis")
            print("y           - Display slice along y axis")
            print("z           - Display slice along z axis")
            print("Z           - Undo last path")
            print("?           - Print this menu")
            print("\n")

    def KeypressCallback(self, caller, ev):
        iren = caller
        key_sym = iren.GetKeySym()
        self.KeypressCommands(key_sym)


# Initialize everything
if __name__ == "__main__":
    # Get the folder and such
    input_paths = GetFileAndFolderNames()
    # input_paths = ['C:/Users/satya/Documents/PycharmProjects/VIA Lab/Test Brain - Copy', 'ExtractedCOW.nii', 'ExtractedCOW']
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

    # Initialize ttGUI class
    ttGUI = TreetrackGUI()

    # Set initial plane orientation and image
    ttGUI.current_plane_orientation = 2
    ttGUI.current_image = 0

    # Create plane widget
    plane_widget = vtk.vtkImagePlaneWidget()
    plane_widget.SetInteractor(render_window_interactor)
    plane_widget.SetTextureInterpolate(False)
    plane_widget.SetResliceInterpolateToNearestNeighbour()
    plane_widget.SetInputData(plane_widget_images[ttGUI.current_image][0])  # Initial Image
    plane_widget.SetPlaneOrientationToZAxes()
    plane_widget.SetMarginSizeX(0.0)
    plane_widget.SetMarginSizeY(0.0)

    # Set initial window and level
    plane_widget.SetWindowLevel(plane_widget_images[ttGUI.current_image][1], plane_widget_images[ttGUI.current_image][2])

    # Set initial plane index
    vtk_dimensions = flipped_vtk_input_image.GetDimensions()
    image_plane_indices = [vtk_dimensions[0] // 2, vtk_dimensions[1] // 2, vtk_dimensions[2] // 2]
    plane_widget.SetSliceIndex(image_plane_indices[ttGUI.current_plane_orientation])

    # Initialize the plane widget
    plane_widget.On()

    # Set up callbacks
    render_window_interactor.AddObserver(vtk.vtkCommand.KeyPressEvent, ttGUI.KeypressCallback)

    # Start the rendering
    render_window.Render()
    render_window_interactor.Initialize()
    render_window_interactor.Start()