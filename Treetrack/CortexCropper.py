from VIA_methods import *

class CortexCropperGUI:
    def __init__(self):
        self.plane_widget_on = True
        self.current_image = None
        self.current_plane_orientation = None

    def KeypressCommands(self, key_sym):

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

        elif key_sym == "x":
            # Store last slice index before switching
            StoreCurrentSliceIndex(plane_widget, image_plane_indices, self.current_plane_orientation, plane_widget_images, self.current_image)
            # Switch to x-axis
            plane_widget.SetPlaneOrientationToXAxes()
            self.current_plane_orientation = 0
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])
            render_window.Render()

        elif key_sym == "D":
            # Toggle on and off the plane widget
            if self.plane_widget_on:
                plane_widget.Off()
            elif not self.plane_widget_on:
                plane_widget.On()

            self.plane_widget_on = not self.plane_widget_on

        elif key_sym == "e":
            # Getting bounds from box widget
            poly_data = vtk.vtkPolyData()
            box_widget.GetPolyData(poly_data)
            bounds = poly_data.GetBounds()

            # Get the min and max bounds
            min_bounds, max_bounds = [0, 0, 0], [0, 0, 0]
            # Convert bounds into ITK-compatible format
            for i in range(0, 6, 2):
                min_bounds[i // 2] = bounds[i]
            for i in range(1, 6, 2):
                max_bounds[i // 2] = bounds[i]

            # Get the min and max indices
            min_index, max_index = ITKIndexAndPhysical(itk_input_image, min_bounds, max_bounds, IndexToPhysical=False)
            # Region size
            region_size = [max_index[i] - min_index[i] for i in range(3)]
            # Create the region for the extract
            region = itk.ImageRegion[3]()
            region.SetIndex(min_index)
            region.SetSize(region_size)

            # Create the ExtractImageFilter
            extract_filter = itk.ExtractImageFilter[itk_input_image, itk_input_image].New()
            extract_filter.SetInput(itk_input_image)
            extract_filter.SetExtractionRegion(region)
            extract_filter.Update()

            # Get the extracted image
            extracted_image = extract_filter.GetOutput()

            # Save the extracted image
            output_path = os.path.join(input_paths[0], input_paths[2] + "_crop.nii")
            itk.imwrite(extracted_image, output_path)

            print("\nCropped and saved successfully as _crop.nii")

        elif key_sym == "m":
            # Getting bounds from box widget
            poly_data = vtk.vtkPolyData()
            box_widget.GetPolyData(poly_data)
            bounds = poly_data.GetBounds()

            # Get the min and max bounds
            min_bounds, max_bounds = [0, 0, 0], [0, 0, 0]
            # Convert bounds into ITK-compatible format
            for i in range(0, 6, 2):
                min_bounds[i // 2] = bounds[i]
            for i in range(1, 6, 2):
                max_bounds[i // 2] = bounds[i]

            # Get the min and max indices
            min_index, max_index = ITKIndexAndPhysical(itk_input_image, min_bounds, max_bounds, IndexToPhysical=False)
            # Region size
            region_size = [max_index[i] - min_index[i] for i in range(3)]
            # Create the region for the extract
            region = itk.ImageRegion[3]()
            region.SetIndex(min_index)
            region.SetSize(region_size)

            # Create the ExtractImageFilter
            extract_filter = itk.ExtractImageFilter[itk_input_image, itk_input_image].New()
            extract_filter.SetInput(itk_input_image)
            extract_filter.SetExtractionRegion(region)
            extract_filter.Update()

            # Get the extracted image
            extracted_image = extract_filter.GetOutput()

            # Save the extracted image
            output_path = os.path.join(input_paths[0], input_paths[2] + "_remap.nii")
            itk.imwrite(extracted_image, output_path)

            print("\nCropped and saved successfully as _remap.nii")

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
            print("i           - Toggle on/off box widget")
            print("D           - Toggle on/off plane widget")
            print("e           - Save current box")
            print("x           - Display slice along x axis")
            print("y           - Display slice along y axis")
            print("z           - Display slice along z axis")
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
    filename = os.path.join(input_paths[0], input_paths[1])
    print(f"Working on brain scan: {input_paths[2]}")

    # Read NIfTI file using ITK
    itk_input_image = itk.imread(filename)

    # Print size
    input_image_size = itk_input_image.GetLargestPossibleRegion().GetSize()
    print(f"Image Dimensions: {input_image_size[0]} x {input_image_size[1]} x {input_image_size[2]}")

    # Reset origin to (0,0,0) and direction matrix to identity
    ResetOriginAndDirection(itk_input_image)

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

    # Initialize ccGUI class
    ccGUI = CortexCropperGUI()

    # Set initial plane orientation and image
    ccGUI.current_plane_orientation = 2
    ccGUI.current_image = 0

    # Create plane widget
    plane_widget = vtk.vtkImagePlaneWidget()
    plane_widget.SetInteractor(render_window_interactor)
    plane_widget.SetTextureInterpolate(False)
    plane_widget.SetResliceInterpolateToNearestNeighbour()
    plane_widget.SetInputData(plane_widget_images[ccGUI.current_image][0])  # Initial Image
    plane_widget.SetPlaneOrientationToZAxes()
    plane_widget.SetMarginSizeX(0.0)
    plane_widget.SetMarginSizeY(0.0)

    # Set initial window and level
    plane_widget.SetWindowLevel(plane_widget_images[ccGUI.current_image][1], plane_widget_images[ccGUI.current_image][2])

    # Get rough COW indices
    itk_min_index, itk_max_index = [0, 0, 0], [input_image_size[0], input_image_size[1], input_image_size[2]]
    x_index = itk_min_index[0] + (itk_max_index[0] - itk_min_index[0]) // 2
    y_index = itk_min_index[1] + (itk_max_index[1] - itk_min_index[1]) // 2
    z_index = itk_min_index[2] + (itk_max_index[2] - itk_min_index[2]) // 2
    # Set initial place widget indices
    image_plane_indices = [x_index, y_index, z_index]
    plane_widget.SetSliceIndex(image_plane_indices[ccGUI.current_plane_orientation])

    # Initialize the plane widget
    plane_widget.On()

    # Create a box widget
    box_widget = vtk.vtkBoxWidget()
    box_widget.SetInteractor(render_window_interactor)
    box_widget.SetPlaceFactor(1)
    box_widget.RotationEnabledOff()

    # Set initial indices for box widget
    vtk_min_coordinate, vtk_max_coordinate = ITKIndexAndPhysical(itk_input_image, itk_min_index, itk_max_index, IndexToPhysical=True)
    box_widget.PlaceWidget(vtk_min_coordinate[0], vtk_max_coordinate[0], vtk_min_coordinate[1], vtk_max_coordinate[1], vtk_min_coordinate[2], vtk_max_coordinate[2])

    # Set up callbacks
    render_window_interactor.AddObserver(vtk.vtkCommand.KeyPressEvent, ccGUI.KeypressCallback)

    # Start the rendering
    render_window.Render()
    render_window_interactor.Initialize()
    render_window_interactor.Start()
