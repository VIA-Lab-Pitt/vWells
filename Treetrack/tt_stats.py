import os
from collections import defaultdict

class vesselResults:
    def __init__(self, FromSaveObject):
        self.FromSave = FromSaveObject

        # Actual results
        self.data_dict = defaultdict(dict)
        self.vessels = ("BA", "L-PCA", "R-PCA", "L-PcomA", "R-PcomA", "L-ICA", "R-ICA", "L-ACA", "R-ACA", "AcomA")

    def AddData(self, vessel_index, key, data):
        current_vessel = self.vessels[vessel_index]
        self.data_dict[current_vessel][key] = data

    def ResetDataFrame(self):
        self.data_dict = defaultdict(dict)

    def FinalBrainCalculations(self):
        pass

    def SaveResults(self, input_paths, resampling_factor, file_extension):
        # New Folder for results stuff
        input_paths.append(f"{input_paths[2]}{file_extension}{resampling_factor}")
        # Check if the new directory exists
        directory_path = os.path.join(input_paths[0], input_paths[-1])
        if not os.path.exists(directory_path):
            # If not, create the directory
            os.makedirs(directory_path)

        # Save the default dict and Excel file
        print(f"\nSaving results to: {directory_path}")
        self.FromSave.SaveDefaultDicts(directory_path, file_extension, [self.data_dict])
        self.FromSave.SaveExcelFile(directory_path, file_extension, self.data_dict)
        print("\tResults have been saved successfully.")
