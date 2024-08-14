import numpy as np
import pandas as pd
from protein import *
from dataclasses import dataclass
import matplotlib.pyplot as plt
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objects as go
import os
from variables import *

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
                with open(f'peaklists/{name}_empty.csv', 'r') as f:
                    print(f'peaklists/{name}_empty.csv loaded')
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

    def get_assigned_peaks(self, fully = False, dim_index: int = None):
        # check if any or all columns starting with Assign is not None
        assignment_columns = self.get_assignment_column_names()
        if fully:
            assigned_peaks = self.peaklist[self.peaklist[assignment_columns].notna().all(axis=1)]
        elif dim_index is not None:
            assigned_peaks = self.peaklist[self.peaklist[assignment_columns[dim_index]].notna()]
        else:
            assigned_peaks = self.peaklist[self.peaklist[assignment_columns].notna().any(axis=1)]
        return PeakList(self.name + '_assigned', self.dimensions, self.spectrum, assigned_peaks)

    def get_unassigned_peaks(self, fully = False, dim_index: int = None):
        assignment_columns = self.get_assignment_column_names()
        if fully:
            unassigned_peaks = self.peaklist[~self.peaklist[assignment_columns].isna().all(axis=1)]
        elif dim_index is not None:
            unassigned_peaks = self.peaklist[~self.peaklist[assignment_columns[dim_index]].isna()]
        else:
            unassigned_peaks = self.peaklist[~self.peaklist[assignment_columns].notna().any(axis=1)]
        return PeakList(self.name + '_unassigned', self.dimensions, self.spectrum, unassigned_peaks)

    def get_peaks_along_dimension(self, location: tuple, dim_index: int, tolerance = 1, num_peaks = None, cs_range: tuple = None):
        """
        Get all peaks along a dimension with the given location and tolerance.
        The resulting peaklist matches the given location in all dimensions except the one specified by dim_index.
        """
        peaklist = self.peaklist
        # filter the peaklist by the given location and tolerance
        if tolerance:
            condition = [] # a list of booleans indicating whether a row satisfies the condition
            for row in range(len(peaklist)):
                if all([abs(peaklist[f'Position F{i+1}'][row] - location[i]) < tolerance for i in range(len(location)) if i != dim_index]):
                    condition.append(True)
                else:
                    condition.append(False)

            peaks = peaklist[condition]
        else:
            peaks = peaklist

        # filter by range if specified
        if cs_range:
            min_val, max_val = cs_range
            peaks = peaks[(peaks[f'Position F{dim_index+1}'] >= min_val) & (peaks[f'Position F{dim_index+1}'] <= max_val)]

        if peaks.empty:
            print('No peaks found.')
            return

        # sort the peaks by lowest MAE
        peaks['MAE'] = peaks.index.map(lambda row: np.mean([abs(peaks[f'Position F{i+1}'][row] - location[i]) for i in range(len(location)) if i != dim_index]))
        sorted_peaks = peaks.sort_values(by='MAE')

        # if num_peaks is specified, return the top num_peaks peaks by lowest MAE
        if num_peaks and num_peaks < len(peaks):
            print(f'Returning {num_peaks} out of {len(peaks)} peaks with lowest MAE:')
            print(sorted_peaks[:num_peaks])
            return sorted_peaks[:num_peaks]

        # otherwise, return sorted peaks
        return sorted_peaks


    def update_assignment(self, location: tuple, dimension: int, protein, residue_index, residue_name, atom):
        """
        Update the assignment of a peak at a given location.
        """
        # dim_index = self.dimensions.index(dimension)
        assignment_str = str(residue_index) + residue_name + atom

        peak = self.get_peak(location)
        print("*********************")
        print(f'peaklists/{self.name}_autoassign_new.csv')

        print(type(self.peaklist[f'Assign F{dimension+1}'][peak.index]))
        print(self.peaklist[f'Assign F{dimension+1}'][peak.index])

        self.peaklist[f'Assign F{dimension+1}'][peak.index] = assignment_str
        print('Peak:\n', peak)

        # self.save_to_csv(f'peaklists/{self.name}_autoassign.csv')

        self.save_to_csv(f'peaklists/{self.name}_autoassign_new.csv')

    def assign_peak(self, peak, dimension: int, protein, residue_index, residue_name, atom):
        """
        Assign an atom to the given peak.
        """
        assignment_str = str(residue_index) + residue_name + atom
        previous_assignment = self.peaklist.at[peak.index[0], f'Assign F{dimension+1}'] # find the previous assignment
        # print("Previous assignment:", previous_assignment)

        # if previous_assignment:
        #     print(previous_assignment)
        #     # print(f'Warning: Overwriting previous assignment {previous_assignment} for peak at {location} in dimension {dimension+1}')
        #     raise ValueError(f'Cannot overwrite previous assignment {previous_assignment}. Exiting...')

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
            text=assigned.peaklist[self.get_assignment_column_names()].apply(lambda row: ', '.join(row.values.astype(str)), axis=1),
            name='Assigned',
            hoverinfo='x+y+z+text'
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
    assert protein[residue_index].three_letter_code == residue_name, f"Residue {residue_index} is not {residue_name}!"
    peaklist.update_assignment(location, dimension, protein, residue_index, residue_name, atom)
    assignment = Assignment(location[dimension], tentative=True)
    protein.assign_atom(residue_index, atom, assignment)

    print(f"\nAssigning to atom: {residue_index} {residue_name}, {atom}...\t", protein[residue_index].get_assignment(atom))


def extend_assignment_left(protein: Protein, residue_index):
    """
    Extend the assignment of a protein to the left by one residue.
    Maps the given N and CA chemical shifts to the CONCA peaklist, then finds the closest aligned peak in CONCA.
    Then the C dimension of the chosen peak is assigned to the i-1 residue.
    """
    # residue_index = residue_index -1
    current_residue = protein[residue_index]
    print("--------------------------------------------------")
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

    # find the closest peak in CONCA peaklist
    peak = CONCA.get_peaks_along_dimension((CA, N), 2, tolerance=TOLERANCE, num_peaks=1, cs_range=target_range)
    if peak is None:
        print('No peaks found in CONCA peaklist matching the target range. Exiting...')
        return

    prev_C = peak['Position F3']
    print(f'Generated assignment for previous residue:', prev_C.values[0])
    assignment = Assignment(prev_C.values[0], tentative=True)
    print(get_location(peak, CONCA.num_dims))
    # assign(CONCA, get_location(peak, CONCA.num_dims), 2, protein, residue_index - 1, prev_residue.three_letter_code, 'C')
    CONCA.assign_peak(peak, 2, protein, residue_index - 1, prev_residue.three_letter_code, 'C')
    prev_residue.assign_atom('C', assignment)

    print(f"\nAssignment extended to {residue_index - 1} {prev_residue.three_letter_code}.")
    print(protein[residue_index - 1])

    # CONCA.plot_peaks_interactive_3d()


def extend_assignment_right(protein: Protein, residue_index):
    pass


def find_N(protein: Protein, residue_index):
    """
    Complete the assignment of an amino acid by assigning its N chemical shift.
    Looks in NCACX for the N value that match all the assigned Cx values.
    """
    print("--------------------------------------------------")
    print(f"Finding N chemical shift for {residue_index} {protein[residue_index].three_letter_code}...")

    C = protein[residue_index].get_assignment('C').chemical_shift
    print('C:', C)
    CA = protein[residue_index].get_assignment('CA').chemical_shift
    print('CA:', CA)

    if C is None or CA is None:
        raise ValueError('Unable to find N chemical shift. Missing assignment for C or CA!')
    NCACX = PeakList('NCACX', ('13C', '13C', '15N'), filename='peaklists/NCACX_autoassign.csv')

    # try N values that match all the assigned Cx values
    min_n, max_n = 105, 135 # assuming N values fall between 105 and 135
    step_size = 0.1

    def match_cx_values(N, peaks, Cx_assignments):
        Cx_guesses = peaks["Position F1"].to_list()
        # loop over current assignments and check if Cx values match the peaks found in the N dimension
        for atom, assignment in Cx_assignments.items():
            if not assignment:
                continue # skip over atom if assignment is None
            # check if Cx value matches any of the peak found in the N dimension
            Cx_real = assignment.chemical_shift
            # peaks = peaks[peaks["Position F1"] in range(N - TOLERANCE, N + TOLERANCE)]
            if all([abs(Cx_real - Cx_guess) > TOLERANCE for Cx_guess in Cx_guesses]):
                print(f"Cx chemical shift {Cx_real} not found in this N dimension.")
                return False # Cx values do not match for this N value
            # print(f"N chemical shift found: {N}")
            # protein[residue_index].assign_atom('N', Assignment(N, tentative=True))
        print(f"N chemical shift found: {N}")
        return True # Cx values match for this N value

    possible_N_values = [] # list to store possible N values
    N = min_n # initial guess

    while N <= max_n:
        print('Trying N:', N)
        peaks = NCACX.get_peaks_along_dimension((C, CA, N), 0, tolerance=TOLERANCE, num_peaks=None, cs_range=None)
        if peaks is not None:
            print('Found possible N chemical shift:', N, "\nChecking Cx values...')")
            if match_cx_values(N, peaks, protein[residue_index].get_assignments()):
                print(peaks)
                # protein[residue_index].assign_atom('N', Assignment(N, tentative=True))
                possible_N_values.append(N)
        N += step_size
    print(possible_N_values)
    print(protein[residue_index])
    if not possible_N_values:
        raise ValueError(f'Unable to find N chemical shift! Please check the assignments for this residue: {residue_index}{protein[residue_index].three_letter_code}{protein[residue_index]}')
    return possible_N_values


def complete_assignment(protein: Protein, residue_index):
    """
    Complete the assignment of an amino acid.
    The C atom and N(i+1) should already be assigned before calling this function.
    The generated assignments should include the Cx and N values.
    """
    print("--------------------------------------------------")

    C = protein[residue_index].get_assignment('C').chemical_shift
    next_N = protein[residue_index + 1].get_assignment('N').chemical_shift

    print(f"Completing assignment for {residue_index} {protein[residue_index].three_letter_code}...")
    print('C:', C, 'N(i+1):', next_N)

    if C is None or next_N is None:
        raise ValueError('Unable to complete assignment. Missing assignments for C or N(i+1)!')

    NCOCX = PeakList('NCOCX', ('13C', '13C', '15N'), filename='peaklists/NCOCX_autoassign.csv')

    target_ranges = AminoAcid.get_cx_chemical_shift_ranges(protein[residue_index].three_letter_code)
    print('Target ranges:', target_ranges)

    peaks = NCOCX.get_peaks_along_dimension((C, C, next_N), 0, tolerance=TOLERANCE, num_peaks=None, cs_range=None)
    print_filtered(peaks)

    for atom, range in target_ranges.items():
        min, max = range
        # pick the peak with the highest intensity within the specified range
        peak = peaks.loc[(peaks["Position F1"] >= range[0]) & (peaks["Position F1"] <= range[1])]
        if peak.empty:
            print(f"No {atom} peak found within range {range}!")
            continue
        print_filtered(peak)
        # peak = peak.sort_values(atom, ascending=False).iloc[0]
        peak = peak.sort_values('Volume', ascending=False).iloc[0]
        print(f"Selected peak for {atom}:")
        print_filtered(peak)
        location = get_location(peak, NCOCX.num_dims, 0)
        assignment = Assignment(location, tentative=True)
        protein[residue_index].assign_atom(atom, assignment)
        print(f"Assigned {atom} to {location}...")

    print("All Cx atoms assigned successfully!")
    print("Finding N in NCACX...")
    # NCACX = PeakList('NCACX', ('13C', '15N', '13C'), filename='peaklists/NCACX_autoassign.csv')
    # target_range = AminoAcid.get_chemical_shift_range(protein[residue_index].three_letter_code, 'N')
    # peaks = NCACX.get_peaks_along_dimension((C, C, next_N), 0, tolerance=1, num_peaks=None, cs_range=target_range)

    possible_N_values = find_N(protein, residue_index)
    for N in possible_N_values:
        try:
            protein[residue_index].assign_atom('N', Assignment(N, tentative=True))
            extend_assignment_left(protein, residue_index)
            complete_assignment(protein, residue_index - 1)
        except ValueError as e:
            print(f"Error assigning N chemical shift: {e}")
            continue

    # print(f"\nAssignment completed for {residue_index} {protein[residue_index].three_letter_code}.")
    # print(protein[residue_index])





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
    try:
        columns = [
            col for col in df.columns if 'Assign' in col or 'Position' in col
        ]
        print(df[columns])
        return df[columns]
    except AttributeError:
        print(df)
        return df


if __name__ == "__main__":
    pass
