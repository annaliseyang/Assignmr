import numpy as np
import pandas as pd
from protein import *
from dataclasses import dataclass
import matplotlib.pyplot as plt
import plotly.express as px

@dataclass
class Assignment:
    chemical_shift: float
    tentative: bool

# initialize peaklists
class PeakList:
    def __init__(self, name: str, dimensions: tuple, spectrum = None, peaklist: pd.DataFrame = None, filename = None) -> None:
        self.name = name
        self.dimensions = dimensions
        self.spectrum = spectrum

        if peaklist is not None:
            self.peaklist = peaklist
        elif filename is not None:
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
        return '\n' + self.name + '\n' + str(df.head())

    def get_position_column_names(self):
        return [f'Position F{i+1}' for i in range(len(self.dimensions))]

    def position_columns(self):
        return self.peaklist[self.get_position_column_names()]

    def get_assignment_column_names(self):
        return [f'Assign F{i+1}' for i in range(len(self.dimensions))]

    def assignment_columns(self):
        return self.peaklist[self.get_assignment_column_names()]

    def get_peak(self, location: tuple):
        return self.peaklist[(self.peaklist['Position F1'] == location[0]) & (self.peaklist['Position F2'] == location[1])]

    def get_assigned_peaks(self, fully = False):
        # check if any or all columns starting with Assign is not None
        assignment_columns = self.get_assignment_column_names()
        if fully:
            assigned_peaks = self.peaklist[self.peaklist[assignment_columns].notna().all(axis=1)]
        else:
            assigned_peaks = self.peaklist[self.peaklist[assignment_columns].notna().any(axis=1)]
        return PeakList(self.name + '_assigned', self.dimensions, self.spectrum, assigned_peaks)

    def get_unassigned_peaks(self, fully = False):
        assignment_columns = self.get_assignment_column_names()
        if fully:
            unassigned_peaks = self.peaklist[~self.peaklist[assignment_columns].isna().all(axis=1)]
        else:
            unassigned_peaks = self.peaklist[~self.peaklist[assignment_columns].notna().any(axis=1)]
        return PeakList(self.name + '_unassigned', self.dimensions, self.spectrum, unassigned_peaks)

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
        """
        Update the assignment of a peak at a given location.
        """
        dim_index = self.dimensions.index(dimension)
        assignment_str = str(residue_index) + residue_name + atom

        peak = self.get_peak(location)
        self.peaklist[f'Assign F{dim_index+1}'][peak.index] = assignment_str
        print('Peak:\n', peak)

        self.save_to_csv(f'peaklists/{self.name}_autoassign.csv')

    def save_to_csv(self, filename):
        self.peaklist.to_csv(filename, index=False)

    def plot_peaks(self, assigned_only = False, annotations = False):
        """
        Plot the peaks in a peaklist.
        """
        assigned = self.get_assigned_peaks()

        fig, ax = plt.subplots(figsize=(10, 7))

        if not assigned_only:
            sc = ax.scatter(self.peaklist['Position F1'], self.peaklist['Position F2'], c=self.peaklist['Volume'], cmap='gray', s=10, alpha=0.4)
        sc = ax.scatter(assigned.peaklist['Position F1'], assigned.peaklist['Position F2'], c=assigned.peaklist['Volume'], cmap='magma')

        plt.colorbar(sc, ax=ax, label='Intensity')

        if annotations:
            for i, row in assigned.assignment_columns().iterrows():
                txt = row.values
                ax.annotate(txt, (self.peaklist['Position F1'][i], self.peaklist['Position F2'][i]))

        ax.set_xlabel(f'Position F1 ({self.dimensions[0]})')
        ax.set_ylabel(f'Position F2 ({self.dimensions[1]})')
        ax.set_title(f'{self.name} Peaks')
        plt.show()

    def plot_peaks_interactive(self, assigned_only = False):
        assigned = self.get_assigned_peaks()
        unassigned = self.get_unassigned_peaks()

        fig = px.scatter(assigned.peaklist, x='Position F1', y='Position F2', color='Volume', size_max=100, opacity=1, color_continuous_scale='magma', hover_data=self.get_assignment_column_names())
        if not assigned_only:
            fig.add_trace(px.scatter(unassigned.peaklist, x='Position F1', y='Position F2', color='Volume', size_max=10, opacity=0.3, color_continuous_scale='greys').data[0])
        fig.show()


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
