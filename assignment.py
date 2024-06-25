import numpy as np
import pandas as pd
from protein import *
from dataclasses import dataclass
import matplotlib.pyplot as plt
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objects as go
import os

@dataclass
class Assignment:
    chemical_shift: float
    tentative: bool

# initialize peaklists
class PeakList:
    def __init__(self, name: str, dimensions: tuple, spectrum = None, peaklist: pd.DataFrame = None, filename = None) -> None:
        self.name = name
        self.dimensions = dimensions
        self.num_dims = len(dimensions)
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
            print(sorted_peaks[:num_peaks])
            return sorted_peaks[:num_peaks]

        # otherwise, return all peaks
        return peaks


    def update_assignment(self, location: tuple, dimension: int, protein, residue_index, residue_name, atom):
        """
        Update the assignment of a peak at a given location.
        """
        # dim_index = self.dimensions.index(dimension)
        assignment_str = str(residue_index) + residue_name + atom

        peak = self.get_peak(location)
        self.peaklist[f'Assign F{dimension+1}'][peak.index] = assignment_str
        print('Peak:\n', peak)

        self.save_to_csv(f'peaklists/{self.name}_autoassign.csv')

    def assign_peak(self, peak, dimension: int, protein, residue_index, residue_name, atom):
        """
        Assign an atom to the given peak.
        """
        assignment_str = str(residue_index) + residue_name + atom
        self.peaklist[f'Assign F{dimension+1}'][peak.index] = assignment_str

        # print the updated peak
        new_peak = self.peaklist.loc[peak.index]

        print(f'New assignment made to {self.name} peak:')
        print_filtered(new_peak)

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
        fig.update_traces(marker_size=10)
        if not assigned_only:
            fig.add_trace(px.scatter(unassigned.peaklist, x='Position F1', y='Position F2', color='Volume', size_max=10, opacity=0.3, color_continuous_scale='greys').data[0])
            fig.update_traces(selector={'opacity': 0.3}, marker_size=5)

        fig.update_xaxes(title_text=f'Position F1 ({self.dimensions[0]})')
        fig.update_yaxes(title_text=f'Position F2 ({self.dimensions[1]})')
        fig.update_layout(title=self.name)

        os.makedirs('plots', exist_ok=True)
        plot(fig, filename=f'plots/{self.name}.html')

    def plot_peaks_interactive_3d(self, assigned_only = False):
        assigned = self.get_assigned_peaks()
        unassigned = self.get_unassigned_peaks()

        fig = go.Figure()
        fig.add_trace(go.Scatter3d(
            x=assigned.peaklist['Position F1'],
            y=assigned.peaklist['Position F2'],
            z=assigned.peaklist['Position F3'],
            mode='markers',
            marker=dict(
                size=5,
                color=assigned.peaklist['Volume'],
                colorscale='magma',
                opacity=1
            ),
            text=self.get_assignment_column_names(),
            name='Assigned'
        ))

        if not assigned_only:
            fig.add_trace(go.Scatter3d(
            x=unassigned.peaklist['Position F1'],
            y=unassigned.peaklist['Position F2'],
            z=unassigned.peaklist['Position F3'],
            mode='markers',
            marker=dict(
                size=2.5,
                color=unassigned.peaklist['Volume'],
                colorscale='gray',
                opacity=0.15
            ),
            text=self.get_assignment_column_names(),
            name='Unassigned'
            ))

        fig.update_layout(
            title={'text': f'{self.name} 3D'},
            scene=dict(
                xaxis_title=f'Position F1 ({self.dimensions[0]})',
                yaxis_title=f'Position F2 ({self.dimensions[1]})',
                zaxis_title=f'Position F3 ({self.dimensions[2]})'
            )
        )

        os.makedirs('plots', exist_ok=True)
        plot(fig, filename=f'plots/{self.name}_3d.html')


def assign(peaklist: PeakList, location: tuple, dimension: int, protein: Protein, residue_index, residue_name, atom):
    peaklist.update_assignment(location, dimension, protein, residue_index, residue_name, atom)
    assignment = Assignment(location[dimension], tentative=True)
    protein.assign_atom(residue_index, atom, assignment)

    print(f"\nAssigning to atom: {residue_index} {residue_name}, {atom}...\t", protein[residue_index].get_assignment(atom))


def extend_assignment_left(protein: Protein, residue_index):
    """
    Extend the assignment of a protein to the left by one residue.
    """
    current_residue = protein[residue_index]
    print(f"Extending assignment to the left of {residue_index} {current_residue.three_letter_code}...")
    print('Current residue:', current_residue)

    N = current_residue.get_assignment('N').chemical_shift
    CA = current_residue.get_assignment('CA').chemical_shift

    print('N:', N, 'CA:', CA)

    if N is None or CA is None:
        raise ValueError('Unable to extend assignment. Missing assignments for N or CA!')

    # get the previous residue
    prev_residue = protein[residue_index - 1]
    print('Previous residue:', prev_residue)

    # find peaks in CONCA peaklist given N and CA chemical shifts of the current residue
    CONCA = PeakList('CONCA', ('13C', '15N', '13C'), filename='peaklists/CONCA_autoassign.csv')

    target_range = AminoAcid.get_chemical_shift_range(prev_residue.three_letter_code, 'C')
    print('Target range:', target_range)
    peak = CONCA.get_peaks_along_dimension((CA, N), 2, tolerance=1, num_peaks=1, cs_range=target_range)

    prev_C = peak['Position F3']
    print(f'Generated assignment for previous residue:', prev_C.values[0])
    assignment = Assignment(prev_C.values[0], tentative=True)
    print(get_location(peak, CONCA.num_dims))
    # assign(CONCA, get_location(peak, CONCA.num_dims), 2, protein, residue_index - 1, prev_residue.three_letter_code, 'C')
    CONCA.assign_peak(peak, 2, protein, residue_index - 1, prev_residue.three_letter_code, 'C')

    print(f"\nAssignment extended to {residue_index - 1} {prev_residue.three_letter_code}.")
    print(protein[residue_index - 1])



def get_location(peak, num_dims, dimension: int = None):
    """
    Get the location (ppm) of a peak along a dimension if the dimension is specified.
    Otherwise, return a tuple of the coordinates.
    """
    if dimension is None:
        return tuple(peak[f'Position F{i+1}'].values[0] for i in range(num_dims))
    else:
        return peak[f'Position F{dimension+1}']

def print_filtered(df):
    """
    Print only the 'Assign' and 'Position' columns of a dataframe.
    """
    columns = [
        col for col in df.columns if 'Assign' in col or 'Position' in col
    ]
    print(df[columns].head())
    return df[columns]


if __name__ == "__main__":
    pass
