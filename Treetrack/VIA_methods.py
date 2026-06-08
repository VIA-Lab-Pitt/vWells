import numpy as np
import itk
import vtk
import os
import pickle
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tkinter import filedialog
from scipy.optimize import minimize
from scipy.integrate import quad
import time
from scipy.ndimage import zoom
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from collections import defaultdict


def ArrangeTwoSets3Elements(set1, set2):
    # Find the common element
    common_elements = list(set1 & set2)
    if len(common_elements) != 1:
        return False

    # Find the unique elements
    unique_set1 = list(set1 - set2)
    unique_set2 = list(set2 - set1)
    if len(unique_set1) != 1 or len(unique_set2) != 1:
        return False

    # Create the list with the desired order
    return [unique_set1[0], common_elements[0], unique_set2[0]]

def ArcLengthFromFittedPolynomial(coefficients, t_bounds):
    # coefficients holds the polynomial coefficients for x(t), y(t), and z(t)
    # t_bounds is the range of t to integrate over

    # Derivatives of the parametric equations
    def dx_dt(t):
        return np.polyval(np.polyder(coefficients[0]), t)

    def dy_dt(t):
        return np.polyval(np.polyder(coefficients[1]), t)

    def dz_dt(t):
        return np.polyval(np.polyder(coefficients[2]), t)

    # The integrand
    def integrand(t):
        return np.sqrt(dx_dt(t) ** 2 + dy_dt(t) ** 2 + dz_dt(t) ** 2)

    # Calculate the arc length
    arc_length, _ = quad(integrand, t_bounds[0], t_bounds[1])

    return arc_length

def CalculatePointToSegmentDistance(point, segment_start, segment_end):
    # Convert points to numpy arrays for vector operations
    point_position = np.array(point)
    segment_start_position = np.array(segment_start)
    segment_end_position = np.array(segment_end)

    # Calculate vectors
    segment_vector = segment_end_position - segment_start_position
    start_to_point_vector = point_position - segment_start_position
    segment_length_squared = np.dot(segment_vector, segment_vector)

    # Project start_to_point_vector onto the segment
    projection_scalar = np.dot(start_to_point_vector, segment_vector) / segment_length_squared

    # Clamp the projection to the range [0, 1] so it stays on the segment
    if isinstance(projection_scalar, np.ndarray):
        projection_scalar = np.clip(projection_scalar, 0, 1)
    else:
        projection_scalar = max(0, min(1, projection_scalar))

    # Closest point(s) on the segment to the point in question
    closest_point = segment_start_position + projection_scalar * segment_vector

    # Distance from the point to the closest point on the segment
    distance = np.linalg.norm(point_position - closest_point)

    return distance

def CheckFragmentInfoReversed(fragments, pixel_locations, path_of_start_nodes):
    fragment_info_reversed = []
    for i, fragment in enumerate(fragments):
        pixels_in_fragment = pixel_locations.get(fragment)
        start_node = path_of_start_nodes[i]
        # Store True if the pixels are stored reversed relative to the start node, else False
        if pixels_in_fragment[-1] == start_node:
            fragment_info_reversed.append(True)
        elif pixels_in_fragment[0] == start_node:
            fragment_info_reversed.append(False)

    return fragment_info_reversed

def ChooseFile():
    return filedialog.askopenfilename(title='Select Input Image')

def ChooseFolder():
    return filedialog.askdirectory(title='Folder with data')

def CreateSphereActor(point, color, radius, list_of_spheres, renderer):
    # Create a sphere at the selected point
    sphere_actor = GetSphereActor(point, color, radius)
    list_of_spheres.append(sphere_actor)
    renderer.AddActor(sphere_actor)

def ConvertFloatListToBinnedDict(list_of_values, list_of_indices, bin_size):
    max_dist = np.max(list_of_values)
    # Bin edges spaced by bin_size, like [0,1), [1,2), and so on
    bin_edges = np.arange(0, max_dist + bin_size, bin_size)

    # Bin each value to its corresponding bin index
    bin_indices = np.digitize(list_of_values, bins=bin_edges)

    # Group indices by bin
    bins = defaultdict(list)
    for upper_bin, index in zip(bin_indices, list_of_indices):
        bins[upper_bin].append(index)

    return bins

def ConvertITKPointToVTK(itk_point, input_image):
    # Input image things
    input_image_spacing = np.array(input_image.GetSpacing())
    input_image_size = input_image.GetLargestPossibleRegion().GetSize()

    # Flip y-axis point
    flipped_itk_point = itk_point.copy()
    flipped_itk_point[1] = input_image_size[1] - itk_point[1] - 1

    # Convert ITK point to VTK
    vtk_point = flipped_itk_point * input_image_spacing

    return vtk_point

def ConvertNumpyPointToVTK(numpy_point, input_image, resampling_factor):
    k, j, i = numpy_point
    itk_point = np.array([i, j, k])

    # Divide by resampling factor
    resampled_itk_point = itk_point / resampling_factor

    # Input image things
    input_image_spacing = np.array(input_image.GetSpacing())
    input_image_size = np.array(input_image.GetLargestPossibleRegion().GetSize())

    # Flip y-axis point
    flipped_itk_point = np.copy(resampled_itk_point)
    flipped_itk_point[1] = input_image_size[1] - resampled_itk_point[1] - 1

    # Multiplying by spacing
    vtk_point = flipped_itk_point * input_image_spacing

    return list(vtk_point)

def ConvertVTKPointToITKNumpyCoordinates(vtk_position, input_image, resampling_factor):
    # Convert to numpy array
    vtk_point = np.array(vtk_position)

    # Divide by spacing and absolute value
    input_image_spacing = np.array(input_image.GetSpacing())
    itk_point = np.abs(vtk_point / input_image_spacing)

    # Flip y-axis point
    input_image_size = np.array(input_image.GetLargestPossibleRegion().GetSize())
    itk_point[1] = input_image_size[1] - itk_point[1] - 1

    # Resample itk point according to the resampling factor
    i, j, k = (int(np.round((coord * resampling_factor))) for coord in itk_point)

    # Switch to convert to numpy
    rounded_numpy_point = (k, j, i)

    return itk_point, vtk_point, rounded_numpy_point

class DataHandling:
    def __init__(self, input_name, resampling_factor):
        # Variables
        self.image_name = input_name
        self.resampling_factor = resampling_factor

    def GetSaveDefaultDictPath(self, directory_path, extension, iteration):
        filename = self.image_name + f"{extension}{self.resampling_factor}_DD{iteration}.pkl"
        return os.path.join(directory_path, filename)

    def GetSaveNumpyPath(self, directory_path, extension):
        filename = self.image_name + f"{extension}{self.resampling_factor}.npz"
        return os.path.join(directory_path, filename)

    def GetSaveXlsxPath(self, directory_path, extension):
        filename = self.image_name + f"{extension}{self.resampling_factor}.xlsx"
        return os.path.join(directory_path, filename)

    def LoadInSavedDefaultDicts(self, directory_path, extension, total_dicts):
        # Load in default dicts
        loaded_dicts = []
        for number in range(total_dicts):
            dict_path = self.GetSaveDefaultDictPath(directory_path, extension, number)
            if os.path.exists(dict_path):
                with open(dict_path, 'rb') as file:
                    loaded_dict = pickle.load(file)
                    loaded_dicts.append(loaded_dict)
            else:
                print(f"\tDictionary does not exist: {dict_path}")
                return None, False

        return loaded_dicts, True

    def LoadInSavedNumpyArrays(self, directory_path, extension):
        # Check to see if numpy image file exists
        numpy_path = self.GetSaveNumpyPath(directory_path, extension)
        if os.path.exists(numpy_path):
            return np.load(numpy_path, allow_pickle=True), True
        else:
            print(f"\tNumpy saved images do not exist: {numpy_path}")
            return None, False

    def SaveNumpyArrays(self, directory_path, extension, objects):
        # Save the Numpy arrays
        numpy_path = self.GetSaveNumpyPath(directory_path, extension)
        np.savez(numpy_path, *objects)

    def SaveDefaultDicts(self, directory_path, extension, objects):
        # Save the default dicts
        for number, default_dict in enumerate(objects):
            dict_path = self.GetSaveDefaultDictPath(directory_path, extension, number)
            with open(dict_path, 'wb') as file:
                pickle.dump(default_dict, file)

    def SaveExcelFile(self, directory_path, extension, dict_of_dicts):
        # Save the Excel file
        excel_filename = self.GetSaveXlsxPath(directory_path, extension)
        with pd.ExcelWriter(excel_filename, engine='xlsxwriter') as writer:
            for key, data in dict_of_dicts.items():
                # Build a DataFrame from the dict, with each list as a column
                df = pd.DataFrame({k: pd.Series(v) for k, v in data.items()})
                # Write the DataFrame as its own sheet
                df.to_excel(writer, sheet_name=key, index=False)

def CurvatureFromFittedPolynomial(coefficients, t_bounds, number_of_points, mm_per_pixel):
    # coefficients holds the polynomial coefficients for x(t), y(t), and z(t)
    # t_bounds is the range of t to evaluate over

    # First derivatives of the parametric equations
    def dx_dt(t):
        return np.polyval(np.polyder(coefficients[0]), t)

    def dy_dt(t):
        return np.polyval(np.polyder(coefficients[1]), t)

    def dz_dt(t):
        return np.polyval(np.polyder(coefficients[2]), t)

    # Second derivatives of the parametric equations
    def d2x_dt2(t):
        return np.polyval(np.polyder(coefficients[0], 2), t)

    def d2y_dt2(t):
        return np.polyval(np.polyder(coefficients[1], 2), t)

    def d2z_dt2(t):
        return np.polyval(np.polyder(coefficients[2], 2), t)

    # Curvature at a given t
    def curvature(t):
        r_prime = np.array([dx_dt(t), dy_dt(t), dz_dt(t)])
        r_double_prime = np.array([d2x_dt2(t), d2y_dt2(t), d2z_dt2(t)])
        numerator = np.linalg.norm(np.cross(r_prime, r_double_prime))
        denominator = np.linalg.norm(r_prime) ** 3
        return numerator / (denominator * mm_per_pixel)

    # Curvature for each t in t_values
    t_values = np.linspace(t_bounds[0], t_bounds[1], number_of_points)
    curvature_values = [curvature(t) for t in t_values]

    return curvature_values

def DefaultDictToNumpyArray(default_dict):
    return np.array(list(default_dict.items()), dtype=object)

def DistanceAndVectorBetweenPoints(point1, point2):
    # Position
    point1_position = np.array(point1)
    point2_position = np.array(point2)
    # Difference
    vector = point2_position - point1_position
    # Calculate distance
    distance = np.sqrt(np.sum(vector ** 2))

    return distance, vector

def DoubleSampleITKImage(loaded_image):
    print("Resampling input image...")
    start = time.perf_counter_ns()

    # Get resampled ITK image spacing
    input_image_spacing = np.array(loaded_image.GetSpacing())
    resampled_image_spacing = input_image_spacing / 2

    # Convert ITK to NumPy array
    input_image = itk.array_from_image(loaded_image).astype(np.float64)
    image_height, image_width, image_depth = input_image.shape

    # Create the new image
    resampled_image_height = (2 * image_height) - 1
    resampled_image_width = (2 * image_width) - 1
    resampled_image_depth = (2 * image_depth) - 1
    resampled_image = np.zeros((resampled_image_height, resampled_image_width, resampled_image_depth), dtype=np.float64)

    # Original pixels keep their values
    resampled_image[0::2, 0::2, 0::2] = input_image

    # Edges are the mean of the 2 adjacent pixels along each axis
    resampled_image[1::2, 0::2, 0::2] = (input_image[:-1, :, :] + input_image[1:, :, :]) / 2
    resampled_image[0::2, 1::2, 0::2] = (input_image[:, :-1, :] + input_image[:, 1:, :]) / 2
    resampled_image[0::2, 0::2, 1::2] = (input_image[:, :, :-1] + input_image[:, :, 1:]) / 2

    # Faces are the mean of the 4 pixels forming a plane
    resampled_image[1::2, 1::2, 0::2] = (input_image[:-1, :-1, :] + input_image[1:, :-1, :] + input_image[:-1, 1:, :] + input_image[1:, 1:, :]) / 4
    resampled_image[0::2, 1::2, 1::2] = (input_image[:, :-1, :-1] + input_image[:, 1:, :-1] + input_image[:, :-1, 1:] + input_image[:, 1:, 1:]) / 4
    resampled_image[1::2, 0::2, 1::2] = (input_image[:-1, :, :-1] + input_image[1:, :, :-1] + input_image[:-1, :, 1:] + input_image[1:, :, 1:]) / 4

    # Cores are the mean of the 8 pixels forming a cube
    resampled_image[1::2, 1::2, 1::2] = (input_image[:-1, :-1, :-1] + input_image[1:, :-1, :-1] + input_image[:-1, 1:, :-1] + input_image[1:, 1:, :-1] + input_image[:-1, :-1, 1:] + input_image[1:, :-1, 1:] + input_image[:-1, 1:, 1:] + input_image[1:, 1:, 1:]) / 8

    end = time.perf_counter_ns()
    print(f'\tDouble sampling took: {(end - start) / 1_000_000_000:.3f} seconds.')
    return resampled_image, resampled_image_spacing

def ExportConvertLinesToTubes(line_actors, radius=0.1, sides=16):
    tube_actors = []
    for line in line_actors:
        # Get the existing polydata
        polydata = line[1].GetMapper().GetInput()

        # Apply tube filter
        tube_filter = vtk.vtkTubeFilter()
        tube_filter.SetInputData(polydata)
        tube_filter.SetRadius(radius)
        tube_filter.SetNumberOfSides(sides)
        tube_filter.Update()

        # Create new actor
        tube_mapper = vtk.vtkPolyDataMapper()
        tube_mapper.SetInputData(tube_filter.GetOutput())

        tube_actor = vtk.vtkActor()
        tube_actor.SetMapper(tube_mapper)
        # Copy color from the original line
        tube_actor.GetProperty().SetColor(line[1].GetProperty().GetColor())

        tube_actors.append(tube_actor)

    return tube_actors

def ExportConvertPlaneWidgetToPolyData(plane_widget):
    # Pull the polydata out of the slice
    plane_geometry = vtk.vtkImageDataGeometryFilter()
    plane_geometry.SetInputData(plane_widget.GetInput())
    plane_geometry.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(plane_geometry.GetOutput())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1, 1, 1)
    return actor


def ExportSceneToGLTF(render_window, segment_actor, line_actors, sphere_actors, plane_widget, filename="exported_scene.glb"):
    output_dir = "3D Models"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_path = os.path.join(output_dir, filename)
    renderer = render_window.GetRenderers().GetFirstRenderer()
    tube_actors = ExportConvertLinesToTubes(line_actors)
    slice_actor = ExportConvertPlaneWidgetToPolyData(plane_widget)
    temp_renderer = vtk.vtkRenderer()

    # Add actors to the temporary renderer
    for actor in tube_actors:
        temp_renderer.AddActor(actor)
    for sphere_actor in sphere_actors:
        temp_renderer.AddActor(sphere_actor)

    # Add the modified segment actor
    temp_renderer.AddActor(segment_actor)

    temp_render_window = vtk.vtkRenderWindow()
    temp_render_window.AddRenderer(temp_renderer)

    # Export to GLTF
    gltf_exporter = vtk.vtkGLTFExporter()
    gltf_exporter.SetRenderWindow(temp_render_window)
    gltf_exporter.SetFileName(file_path)
    gltf_exporter.Write()

    print(f"Exported scene to {file_path}")

def Fit3DPolynomial(actual_points, max_degree, total_fitted_points, r2_threshold, t_bounds=(0, 1)):
    points = np.array(actual_points)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    def poly_func(coefficients, t):
        poly_x = np.polyval(coefficients[0], t)
        poly_y = np.polyval(coefficients[1], t)
        poly_z = np.polyval(coefficients[2], t)
        return np.vstack((poly_x, poly_y, poly_z)).T

    def residuals(coefficients, t, points):
        return np.sum((poly_func(coefficients, t) - points) ** 2)

    def r2_score(actual, predicted):
        ss_res = np.sum((actual - predicted) ** 2)
        ss_tot = np.sum((actual - np.mean(actual, axis=0)) ** 2)
        if ss_tot == 0:
            return 1.0 if ss_res == 0 else -np.inf
        return 1 - (ss_res / ss_tot)

    # Parameter t represents the position along the polynomial
    t = np.linspace(t_bounds[0], t_bounds[1], len(points))

    best_fit_coefficients = None
    best_r2 = -np.inf
    best_degree = None

    for degree in range(1, max_degree + 1):
        # Not enough points to fit the current degree polynomial
        if len(points) < degree + 1:
            break

        # Initial guess for the polynomial coefficients
        coefficients_x = np.polyfit(t, x, degree)
        coefficients_y = np.polyfit(t, y, degree)
        coefficients_z = np.polyfit(t, z, degree)
        initial_guess = [coefficients_x, coefficients_y, coefficients_z]

        # Boundary conditions, the first and last points must be exact
        def boundary_conditions(coefficients):
            start_diff = poly_func(coefficients, np.array([0]))[0] - points[0]
            end_diff = poly_func(coefficients, np.array([1]))[0] - points[-1]
            return np.sum(start_diff ** 2) + np.sum(end_diff ** 2)

        def objective_function(coefficients_flat):
            coefficients = [
                coefficients_flat[:degree + 1],
                coefficients_flat[degree + 1: 2 * (degree + 1)],
                coefficients_flat[2 * (degree + 1):]
            ]
            return residuals(coefficients, t, points) + 1e6 * boundary_conditions(coefficients)

        coefficients_flat_initial = np.concatenate(initial_guess)
        result = minimize(objective_function, coefficients_flat_initial)

        fitted_coefficients_flat = result.x
        fitted_coefficients = [
            fitted_coefficients_flat[:degree + 1],
            fitted_coefficients_flat[degree + 1: 2 * (degree + 1)],
            fitted_coefficients_flat[2 * (degree + 1):]
        ]

        # Generate fitted points
        fitted_points = poly_func(fitted_coefficients, t)

        # Calculate R^2 score
        r2 = r2_score(points, fitted_points)
        if r2 > best_r2:
            best_r2 = r2
            best_fit_coefficients = fitted_coefficients
            best_degree = degree

        if r2 >= r2_threshold:
            break

    # Generate new points along the fitted polynomial
    t_new = np.linspace(t_bounds[0], t_bounds[1], total_fitted_points)
    new_points = poly_func(best_fit_coefficients, t_new)
    # print(f"\nBest degree: {best_degree} with R^2: {best_r2}")
    # print("\tPolynomial equation for x(t):", PrintPolynomialEquations(best_fit_coefficients[0], 't'))
    # print("\tPolynomial equation for x(t):", PrintPolynomialEquations(best_fit_coefficients[0], 't'))
    # print("\tPolynomial equation for x(t):", PrintPolynomialEquations(best_fit_coefficients[0], 't'))

    return new_points, best_fit_coefficients, t_bounds

def Get3DVessel(path_image, image_spacing, opacity, color=(1.0, 1.0, 1.0)):
    # Convert back to ITK
    itk_shortest_path = itk.image_from_array(path_image)
    # Match spacing to input image
    itk_shortest_path.SetSpacing(image_spacing)
    # Convert to VTK and flip along the y-axis
    vtk_shortest_path = ITKtoVTKandFlipY(itk_shortest_path)

    # Create the marching cubes surface from the path image
    marching_cubes_surface = vtk.vtkMarchingCubes()
    marching_cubes_surface.SetInputData(vtk_shortest_path)
    marching_cubes_surface.SetValue(0, 1)
    marching_cubes_surface.Update()

    # Map the surface
    poly_data_mapper = vtk.vtkPolyDataMapper()
    poly_data_mapper.SetInputConnection(marching_cubes_surface.GetOutputPort())
    # Ignore scalar values and use a solid color
    poly_data_mapper.ScalarVisibilityOff()

    # Build the actor
    segmented_vessel_actor = vtk.vtkActor()
    segmented_vessel_actor.SetMapper(poly_data_mapper)
    segmented_vessel_actor.GetProperty().SetOpacity(opacity)
    segmented_vessel_actor.GetProperty().SetColor(color)

    return segmented_vessel_actor

def GetLineActors(point_locations_dict, resampling_factor, itk_input_image, line_width):
    print("Getting line actors for fragments...")

    # Add all actors to this list
    line_actors = []
    for key in point_locations_dict:
        line_actor = GetLineActorFromNumpyPoints(itk_input_image, point_locations_dict.get(key).copy(), resampling_factor, GetVTKColor(key), line_width)

        line_actors.append([key, line_actor])

    return line_actors

def GetFileAndFolderNames():
    absolute_path = ChooseFile()
    input_filename = os.path.basename(absolute_path)
    input_name, _ = os.path.splitext(input_filename)
    folder_path = os.path.dirname(absolute_path)
    return [folder_path, input_filename, input_name]

def GetLineActorFromNumpyPoints(input_image, numpy_points, resampling_factor, color, line_width):
    vtk_points = vtk.vtkPoints()
    total_points = len(numpy_points)

    # Convert all the points to vtk
    for point in numpy_points:
        vtk_point = ConvertNumpyPointToVTK(point, input_image, resampling_factor)
        vtk_points.InsertNextPoint(vtk_point)

    # Create a polyline to connect the points
    polyLine = vtk.vtkPolyLine()
    polyLine.GetPointIds().SetNumberOfIds(total_points)
    for i in range(total_points):
        polyLine.GetPointIds().SetId(i, i)

    # Create a cell array to store the lines in and add the polyline to it
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyLine)

    # Create a polydata to store everything in
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(vtk_points)
    polyData.SetLines(cells)

    # Setup mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    line_actor = vtk.vtkActor()
    line_actor.SetMapper(mapper)
    line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
    line_actor.GetProperty().SetLineWidth(line_width)

    return line_actor

def GetMousePositionAsITKVTKNumpyCoordinates(interactor, renderer, input_image, resampling_factor):
    picker = vtk.vtkCellPicker()
    x, y = interactor.GetEventPosition()
    picker.Pick(x, y, 0, renderer)
    position = picker.GetPickPosition()

    return ConvertVTKPointToITKNumpyCoordinates(position, input_image, resampling_factor)

def GetNearestActorToMouseAsITKVTKNumpyCoordinates(interactor, renderer, input_image, resampling_factor):
    picker = vtk.vtkPointPicker()
    x, y = interactor.GetEventPosition()
    picker.Pick(x, y, 0, renderer)
    position = picker.GetPickPosition()

    return ConvertVTKPointToITKNumpyCoordinates(position, input_image, resampling_factor)

def GetPathImage(index_image, vWells_in_path, dict_of_vWells_and_pixels):
    # An image to be processed by marching cubes
    image_height, image_width, image_depth = index_image.shape
    path_image = np.zeros([image_height, image_width, image_depth])

    # Go through entire list of paths and get a path image
    list_of_vWells = tuple(vWells_in_path)
    total_vWells_in_path = len(list_of_vWells)
    for i in range(total_vWells_in_path):
        vWell = list_of_vWells[i]
        for x, y, z in dict_of_vWells_and_pixels[vWell]:
            path_image[x, y, z] = 2

    print(f"\tTotal vWells in path: {total_vWells_in_path}")

    return path_image

def GetSphereActor(point, color, radius):
    # Create a sphere
    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetCenter(point[0], point[1], point[2])
    sphereSource.SetRadius(radius)
    sphereSource.SetPhiResolution(100)
    sphereSource.SetThetaResolution(100)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphereSource.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))

    return actor

def GetVTKColor(index):
    colors = [
        'Orange', 'Yellow', 'Green', 'Tomato', 'Gold', 'Lime', 'Cyan', 'SkyBlue',
        'Magenta', 'Salmon', 'Coral', 'Chartreuse', 'Turquoise', 'Orchid', 'HotPink',
        'LightGreen', 'Aquamarine', 'Plum', 'DeepPink', 'PaleTurquoise',
        'Teal', 'ForestGreen', 'GoldenRod', 'SandyBrown', 'Peru',
        'RosyBrown', 'MediumPurple', 'LightSalmon', 'MediumSpringGreen', 'Thistle',
        'Olive', 'Moccasin', 'LightCoral', 'PaleGreen', 'PowderBlue', 'LightCyan',
        'Lavender', 'MediumOrchid', 'LightSeaGreen', 'Wheat', 'LawnGreen',
        'PaleGoldenRod', 'MediumAquamarine', 'LightSkyBlue', 'Bisque', 'LightPink',
        'PaleVioletRed', 'Beige', 'BurlyWood', 'DarkSeaGreen', 'AliceBlue',
        'Gainsboro', 'Azure', 'Khaki', 'LightBlue', 'LightSteelBlue',
        'MediumTurquoise', 'MediumSeaGreen', 'Tan', 'Pink',
        'LightSlateGray', 'MediumSpringGreen', 'DarkKhaki', 'MediumTurquoise',
        'IndianRed', 'MediumSeaGreen', 'Goldenrod', 'OliveDrab', 'DarkOrange',
        'LimeGreen', 'SpringGreen', 'DarkSalmon', 'LightSteelBlue',
        'MediumAquamarine', 'MediumVioletRed', 'DarkGoldenrod',
        'SeaGreen', 'MediumOrchid', 'Tomato',
        'Coral', 'DarkOliveGreen', 'Peru', 'LightGreen',
        'MediumPurple', 'MediumSpringGreen', 'LawnGreen', 'SandyBrown',
        'DarkSeaGreen', 'PaleGreen', 'Chartreuse', 'SeaGreen',
        'Turquoise', 'Violet', 'DarkKhaki']

    return colors[index]

def GetVwellFromPoint(numpy_point, index_image):
    # Find the starting and ending vWells (nodes)
    i, j, k = (numpy_point[0], numpy_point[1], numpy_point[2])
    vWell_ID = int(index_image[i, j, k])

    return vWell_ID

def IncrementImageListIndex(image_list, current_index, up=False, down=False):
    if up:
        if current_index < len(image_list) - 1:
            current_index += 1

    elif down:
        if current_index != 0:
            current_index -= 1

    return current_index

def ITKIndexAndPhysical(itk_image, min_index, max_index, IndexToPhysical):
    # Input image size and spacing
    input_image_spacing = np.array(itk_image.GetSpacing())
    input_image_size = itk_image.GetLargestPossibleRegion().GetSize()
    if IndexToPhysical:
        # Flip y coordinate
        temp_max = max_index[1]
        max_index[1] = input_image_size[1] - min_index[1] - 1
        min_index[1] = input_image_size[1] - temp_max - 1
        # Account for spacing
        min_index = (min_index * input_image_spacing).astype(int).tolist()
        max_index = (max_index * input_image_spacing).astype(int).tolist()
    else:
        # Account for spacing
        min_index = np.round(min_index / input_image_spacing).astype(int).tolist()
        max_index = np.round(max_index / input_image_spacing).astype(int).tolist()
        # Flip y coordinate
        temp_max = max_index[1]
        max_index[1] = input_image_size[1] - min_index[1] - 1
        min_index[1] = input_image_size[1] - temp_max - 1

    return min_index, max_index


def ITKtoNumpyAndResamplePoints(itk_points, resampling_factor):
    resampled_numpy_points = []
    for point in itk_points:
        # Find the starting and ending vWells (nodes)
        k, j, i = (point[0] * resampling_factor, point[1] * resampling_factor, point[2] * resampling_factor)
        resampled_numpy_points.append((i, j, k))

    return resampled_numpy_points

def ITKtoVTKandFlipY(itk_image):
    # Convert to VTK
    pre_flip_vtk_image = itk.vtk_image_from_image(itk_image)

    # Flip the image along the y-axis
    flip_filter = vtk.vtkImageFlip()
    flip_filter.SetFilteredAxis(1)
    flip_filter.SetInputData(pre_flip_vtk_image)
    flip_filter.Update()
    vtk_image = flip_filter.GetOutput()

    return vtk_image

def NumpyToITKandSetSpacing(numpy_image, input_image_spacing):
    # Convert to ITK
    itk_spaced_image = itk.image_from_array(numpy_image)
    # Set Spacing
    itk_spaced_image.SetSpacing(input_image_spacing)

    return itk_spaced_image

def PickCurrentActor(interactor, renderer):
    picker = vtk.vtkPropPicker()
    x, y = interactor.GetEventPosition()
    picker.Pick(x, y, 0, renderer)
    # Check if a prop was picked
    if picker.GetActor():
        return picker.GetActor()
    else:
        print("No actor picked")

def PlotEnumeratedValues(y_values, plot_title):
    # The x-values are just the indices of the y-values
    x_values = list(range(len(y_values)))

    # Plot the graph
    plt.plot(x_values, y_values, marker='o')

    # Label the axes
    plt.xlabel('Indices')
    plt.ylabel('Y values')

    # Add a title
    plt.title(plot_title)
    # Show the plot
    plt.show()

def PrintPolynomialEquations(coefficients, var_name):
    # ex) print("\tPolynomial equation for x(t):", display_polynomial(fitted_coefficients[0], 't'))
    degree = len(coefficients) - 1
    terms = [f"{coeff:.2f}{var_name}^{degree - i}" for i, coeff in enumerate(coefficients)]
    equation = " + ".join(terms)
    equation = equation.replace(f"{var_name}^0", "").replace(f"{var_name}^1", var_name).replace("+ -", "- ")
    return equation

def RemoveLastPointResetActor(picked_nodes, node_points, actors, color):
    last_point = picked_nodes.pop()
    index = node_points.index(last_point)
    last_sphere_actor = actors[index]
    last_sphere_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))

def RemoveLastLineResetActor(picked_fragments, line_actors):
    last_fragment = picked_fragments.pop()
    index = next((index for index, (ID, actor) in enumerate(line_actors) if ID == last_fragment), None)
    last_line_actor = line_actors[index][1]
    last_line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(GetVTKColor(last_fragment)))

def ReplaceVesselActor(renderer, segmented_vessel_actor, have_vessel_actor, vWellObject, vessel_opacity):
    if have_vessel_actor:
        # Remove existing actor
        renderer.RemoveActor(segmented_vessel_actor)
    # Get segmented vessel actor
    segmented_vessel_actor = vWellObject.GetSegmentedVesselActor(vessel_opacity)
    # Add actor to renderer
    renderer.AddActor(segmented_vessel_actor)

    return segmented_vessel_actor

def ResampleITKImage(input_image, resampling_factor):
    print(f"Using outdated resampling for input image...")
    # Turn itk image into a numpy array
    numpy_input_image = itk.array_from_image(input_image).astype(np.uint32)
    image_height, image_width, image_depth = numpy_input_image.shape
    # print(f"\tNumpy Dimensions: {image_height} x {image_width} x {image_depth}")

    # Get resampled ITK image spacing
    input_image_spacing = np.array(input_image.GetSpacing())
    resampled_input_image_spacing = input_image_spacing / resampling_factor

    # Double sample the input image by linear interpolation
    resampled_numpy_image = zoom(input_image, resampling_factor, order=1)

    return resampled_numpy_image, resampled_input_image_spacing

def ResetOriginAndDirection(itk_image):
    # Set the direction matrix to the identity
    identityDirectionMatrix = itk.Matrix[itk.D, 3, 3]()
    identityDirectionMatrix.SetIdentity()
    itk_image.SetDirection(identityDirectionMatrix)

    # Set the origin of the input image to (0,0,0)
    originNew = itk.Point[itk.D, 3]()
    originNew.Fill(0.0)
    itk_image.SetOrigin(originNew)

def SetRememberedSliceIndex(plane_widget, current_index, resampling_factor):
    # Set the new slice index from stored values * resampling factor
    new_slice_index = current_index * resampling_factor
    plane_widget.SetSliceIndex(new_slice_index)

def ShortestDistanceToPath(point_in_question, shortest_path):
    distances = []
    # Go through each set of points in the path and get the shortest distance to the segment
    for i in range(len(shortest_path) - 1):
        point1 = shortest_path[i]
        point2 = shortest_path[i + 1]
        # Get the distance to the segment
        distances.append(CalculatePointToSegmentDistance(point_in_question, point1, point2))

    # Get the shortest distance to the shortest path
    shortest_distance = min(distances)

    return shortest_distance

def Show3DPlot(coords, values):
    # Create a 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the skeleton points with color representing the radius
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=values,
                    cmap='viridis')

    # Add a color bar to show the radii values
    plt.colorbar(sc, label='Radius')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Skeleton with Radii')

    # Update the plot's scale on scroll
    def on_scroll(event):
        scale_factor = 1.1 if event.button == 'up' else 0.9
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_factor, scale_factor, scale_factor, 1]))
        fig.canvas.draw_idle()

    # Connect the scroll event to the on_scroll function
    fig.canvas.mpl_connect('scroll_event', on_scroll)

    # Set aspect ratio so the graph responds to window resizing
    ax.set_box_aspect([1, 1, 1])

    plt.show()

def StoreCurrentSliceIndex(plane_widget, image_plane_indices, current_plane_orientation, plane_widget_images, current_image):
    # Get current index
    current_slice_index = plane_widget.GetSliceIndex()
    # Store the current index divided by the resampling factor
    image_plane_indices[current_plane_orientation] = current_slice_index // plane_widget_images[current_image][3]

def UpdatePlaneWidgetImageWindowLevel(plane_widget, plane_widget_images, new_index):
    plane_widget.SetInputData(plane_widget_images[new_index][0])
    plane_widget.SetWindowLevel(plane_widget_images[new_index][1],
                                plane_widget_images[new_index][2])

def Visualize3DGrid(neighborhood):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    neighborhood = np.array(neighborhood)
    # Indices where values are greater than 0
    filled_positions = np.argwhere(neighborhood > 0)

    for pos in filled_positions:
        x, y, z = pos
        # ax.bar3d(x, z, y, 1, 1, 1, shade=True)  # y and z swapped
        # y and z swapped
        ax.bar3d(x, y, z, 1, 1, 1, shade=True)


    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.set_zlabel('Y')

    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_zticks(range(3))
    ax.set_xlim([0, 3])
    ax.set_ylim([0, 3])
    ax.set_zlim([0, 3])

    plt.show()