import numpy as np
import pandas as pd
from protein import *
from dataclasses import dataclass

@dataclass
class Assignment:
    chemical_shift: float
    tentative: bool

# initialize peaklists
class PeakList:
    def __init__(self, name: str, dimensions: tuple, spectrum = None, filename = None) -> None:
        self.name = name
        self.dimensions = dimensions
        self.spectrum = spectrum

        if filename:
            with open(filename, 'r') as f:
                self.peaklist = pd.read_csv(f)
        else:
            try:
                with open(f'peaklists/{name}.csv', 'r') as f:
                    self.peaklist = pd.read_csv(f)
            except FileNotFoundError:
                raise FileNotFoundError(f'Peaklist {name} not found')


    def __str__(self):
        columns = self.get_position_column_names() + self.get_assignment_column_names()
        df = self.peaklist[columns]
        return str(df.head())

    def get_position_column_names(self):
        return [f'Position F{i+1}' for i in range(len(self.dimensions))]

    def get_assignment_column_names(self):
        return [f'Assign F{i+1}' for i in range(len(self.dimensions))]

    def get_peak(self, location: tuple):
        return self.peaklist[(self.peaklist['Position F1'] == location[0]) & (self.peaklist['Position F2'] == location[1])]

    def get_assigned_peaks(self, fully = False):
        # check if any or all columns starting with Assign is not None
        assignment_columns = self.get_assignment_column_names()
        if fully:
            return self.peaklist[self.peaklist[assignment_columns].notna().all(axis=1)]
        else:
            return self.peaklist[self.peaklist[assignment_columns].notna().any(axis=1)]

    def get_unassigned_peaks(self):
        return self.peaklist[~self.peaklist[self.get_assignment_column_names()].notna().any(axis=1)]

    def get_peaks_along_dimension(self, location: tuple, dim_index: int, tolerance = 0.1, num_peaks = None, cs_range: tuple = None):
        """
        Get all peaks along a dimension with the given location and tolerance.
        The resulting peaklist matches the given location in all dimensions except the one specified by dim_index.
        """
        # filter the peaklist by the given location and tolerance
        condition = []
        for row in range(len(self.peaklist)):
            if all([abs(self.peaklist[f'Position F{i+1}'][row] - location[i]) < tolerance for i in range(len(location)) if i != dim_index]):
                condition.append(True)
            else:
                condition.append(False)

        peaks = self.peaklist[condition]

        # filter by range if specified
        if cs_range:
            min_val, max_val = cs_range
            peaks = peaks[(peaks[f'Position F{dim_index+1}'] >= min_val) & (peaks[f'Position F{dim_index+1}'] <= max_val)]

        # if num_peaks is specified, return the top num_peaks peaks by lowest MAE
        if num_peaks and num_peaks < len(peaks):
            peaks['MAE'] = peaks.index.map(lambda row: np.mean([abs(peaks[f'Position F{i+1}'][row] - location[i]) for i in range(len(location)) if i != dim_index]))
            sorted_peaks = peaks.sort_values(by='MAE')
            print(f'Returning {num_peaks} out of {len(peaks)} peaks with lowest MAE:')
            return sorted_peaks[:num_peaks]

        # otherwise, return all peaks
        return peaks


    def update_assignment(self, location: tuple, dimension: str, protein, residue_index, residue_name, atom):
        dim_index = self.dimensions.index(dimension)
        assignment_str = str(residue_index) + residue_name + atom

        peak = self.get_peak(location)
        self.peaklist[f'Assign F{dim_index+1}'][peak.index] = assignment_str
        print('Peak:\n', peak)

        self.save_to_csv(f'peaklists/{self.name}_autoassign.csv')

    def save_to_csv(self, filename):
        self.peaklist.to_csv(filename, index=False)


def assign(peaklist: PeakList, location: tuple, dimension: str, protein: Protein, residue_index, residue_name, atom):
    peaklist.update_assignment(location, dimension, protein, residue_index, residue_name, atom)
    assignment = Assignment(location, True)
    protein.assign_atom(residue_index, atom, assignment)


def print_filtered(df):
    """
    Print only the 'Assign' and 'Position' columns of a dataframe.
    """
    columns = [
        col for col in df.columns if 'Assign' in col or 'Position' in col
    ]
    print(df[columns].head())


if __name__ == "__main__":
    pass
