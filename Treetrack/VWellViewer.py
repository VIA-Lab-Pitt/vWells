import os
from tt_algorithm import *
from VIA_methods import *

class TreetrackGUI:
    def __init__(self):
        self.current_plane_orientation = None
        self.plane_widget_on = True
        self.have_vessel_actor = False
        self.displaying_vessel_actor = False
        self.segmented_vessel_actor = None
        self.vessel_opacity = 0.5
        self.showing_commands = True

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

        elif key_sym == "D":
            # Toggle on and off the plane widget
            if self.plane_widget_on:
                plane_widget.Off()
            elif not self.plane_widget_on:
                plane_widget.On()

            self.plane_widget_on = not self.plane_widget_on

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
            # Switch to z-axis
            plane_widget.SetPlaneOrientationToZAxes()
            self.current_plane_orientation = 2
            # Set the new slice index from stored values
            SetRememberedSliceIndex(plane_widget, image_plane_indices[self.current_plane_orientation], plane_widget_images[self.current_image][3])
            render_window.Render()

        elif key_sym == "question":
            # Toggle the on-screen command lists
            self.showing_commands = not self.showing_commands
            keybind_text_actor.SetVisibility(self.showing_commands)
            mouse_text_actor.SetVisibility(self.showing_commands)
            render_window.Render()

    def KeypressCallback(self, caller, ev):
        iren = caller
        key_sym = iren.GetKeySym()
        self.KeypressCommands(key_sym)


# Initialize everything
if __name__ == "__main__":
    # SETTINGS
    resampling_factor = 2
    kernel_radius = 1

    # Get the folder and such
    input_paths = GetFileAndFolderNames()

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

    # Load/compute vWells (builds Segmentation object, returns variance + mean images)
    itk_variance_image, max_variance, itk_mean_image = vWell_algorithm.StartVwellFindingAlgorithm(input_paths, itk_input_image, resampling_factor, kernel_radius, build_segment_object=True)

    # Mean image
    flipped_vtk_mean_image = ITKtoVTKandFlipY(itk_mean_image)
    plane_widget_images.append([flipped_vtk_mean_image, 230, 95, resampling_factor])

    # Variance image
    flipped_vtk_variance_image = ITKtoVTKandFlipY(itk_variance_image)
    plane_widget_images.append([flipped_vtk_variance_image, 0.09 * max_variance, 0.05 * max_variance, resampling_factor])

    # Load saved segmentation
    if vWell_algorithm.Segmentation.SaveLoadProgress(load=True):
        ttGUI.segmented_vessel_actor = ReplaceVesselActor(renderer, ttGUI.segmented_vessel_actor, ttGUI.have_vessel_actor, vWell_algorithm.Segmentation, ttGUI.vessel_opacity)
        ttGUI.have_vessel_actor = True
        ttGUI.displaying_vessel_actor = True
        print("Loaded existing segmentation.")
    else:
        print("No saved segmentation found.")

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

    # Keybind and mouse command overlay
    min_font_size = 9
    font_fraction = 0.02  # font size as a fraction of window height

    # Keybinds (top-left)
    keybind_text = (
        "Up / Down    - Slice up / down\n"
        "Left / Right - Switch image (input / vWell mean / variance)\n"
        "x / y / z    - View along axis\n"
        "D            - Plane widget on/off\n"
        "V            - Vessel on/off\n"
        "o            - Vessel opacity\n"
        "?            - Toggle this list"
    )
    keybind_text_actor = vtk.vtkTextActor()
    keybind_text_actor.SetInput(keybind_text)
    ktprop = keybind_text_actor.GetTextProperty()
    ktprop.SetFontFamilyToCourier()
    ktprop.SetColor(0.3, 0.6, 1.0)
    ktprop.ShadowOn()
    ktprop.SetJustificationToLeft()
    ktprop.SetVerticalJustificationToTop()
    keybind_text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    keybind_text_actor.SetPosition(0.01, 0.99)
    renderer.AddViewProp(keybind_text_actor)

    # Mouse controls (bottom-left)
    mouse_text = (
        "L-drag (off plane) - Rotate\n"
        "M-drag (off plane) - Translate\n"
        "R-drag (off plane) - Zoom\n"
        "M-drag (on plane) - Slide plane\n"
        "R-drag (on plane) - Brightness/contrast\n"
        "Wheel - Zoom"
    )
    mouse_text_actor = vtk.vtkTextActor()
    mouse_text_actor.SetInput(mouse_text)
    mtprop = mouse_text_actor.GetTextProperty()
    mtprop.SetFontFamilyToCourier()
    mtprop.SetColor(0.3, 0.6, 1.0)
    mtprop.ShadowOn()
    mtprop.SetJustificationToLeft()
    mtprop.SetVerticalJustificationToBottom()
    mouse_text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    mouse_text_actor.SetPosition(0.01, 0.01)
    renderer.AddViewProp(mouse_text_actor)

    # Scale both overlays with window height, clamped to a minimum
    def ScaleCommandFont(caller=None, ev=None):
        window_height = render_window.GetSize()[1]
        font_size = max(min_font_size, round(window_height * font_fraction))
        ktprop.SetFontSize(font_size)
        mtprop.SetFontSize(font_size)
        if caller is not None:
            render_window.Render()

    ScaleCommandFont()
    render_window_interactor.AddObserver(vtk.vtkCommand.ConfigureEvent, ScaleCommandFont)

    # Set up callbacks
    render_window_interactor.AddObserver(vtk.vtkCommand.KeyPressEvent, ttGUI.KeypressCallback)

    # Start the rendering
    render_window.Render()
    render_window_interactor.Initialize()
    render_window_interactor.Start()