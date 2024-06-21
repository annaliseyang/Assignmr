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
        return str(self.peaklist.head())

    def get_position_column_names(self):
        return [f'Position F{i+1}' for i in range(len(self.dimensions))]

    def get_assignment_column_names(self):
        return [f'Assign F{i+1}' for i in range(len(self.dimensions))]

    def get_assigned_peaks(self, fully = False):
        # check if any or all columns starting with Assign is not None
        assignment_columns = self.get_assignment_column_names()
        if fully:
            return self.peaklist[self.peaklist[assignment_columns].notna().all(axis=1)]
        else:
            return self.peaklist[self.peaklist[assignment_columns].notna().any(axis=1)]

    def update_assignment(self, location: tuple, dimension: str, protein, residue_index, residue_name, atom):
        dim_index = self.dimensions.index(dimension)
        assignment_str = str(residue_index) + residue_name + atom

        peak = self.peaklist[(self.peaklist['Position F1'] == location[0]) & (self.peaklist['Position F2'] == location[1])]
        self.peaklist[f'Assign F{dim_index+1}'][peak.index] = assignment_str
        print('Peak:\n', peak)

        self.save_to_csv(f'peaklists/{self.name}_autoassign.csv')

    def save_to_csv(self, filename):
        # filename = f'peaklists/{self.name}_autoassign.csv'
        self.peaklist.to_csv(filename, index=False)


def assign(peaklist: PeakList, location: tuple, dimension: str, protein: Protein, residue_index, residue_name, atom):
    peaklist.update_assignment(location, dimension, protein, residue_index, residue_name, atom)
    protein.assign(residue_index, atom)
