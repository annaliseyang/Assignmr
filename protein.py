import pandas as pd

class AminoAcid:
    __amino_acid_dict = {
        'A': 'Ala',
        'R': 'Arg',
        'N': 'Asn',
        'D': 'Asp',
        'C': 'Cys',
        'Q': 'Gln',
        'E': 'Glu',
        'G': 'Gly',
        'H': 'His',
        'I': 'Ile',
        'L': 'Leu',
        'K': 'Lys',
        'M': 'Met',
        'F': 'Phe',
        'P': 'Pro',
        'S': 'Ser',
        'T': 'Thr',
        'W': 'Trp',
        'Y': 'Tyr',
        'V': 'Val',
    }

    __amino_acid_dict.update( {v: k for k, v in __amino_acid_dict.items()} )

    with open('chemical_shifts.csv', 'r') as f:
        __chemical_shifts = pd.read_csv(f)

    def get_chemical_shifts():
        return AminoAcid.__chemical_shifts

    def __init__(self, id: str):
        id_2 = self.__amino_acid_dict[id]
        match(len(id)):
            case 1:
                self.__one_letter_code = id
                self.__three_letter_code = id_2
            case 3:
                self.__one_letter_code = id_2
                self.__three_letter_code = id

        df_atoms = AminoAcid.__chemical_shifts[ AminoAcid.__chemical_shifts['comp_id'] == self.__three_letter_code ]
        self.__atoms_assignment_states = { atom: False for atom in df_atoms['atom_id'].to_list() }


    def __str__(self) -> str:
        df_atoms = AminoAcid.__chemical_shifts[AminoAcid.__chemical_shifts['comp_id'] == self.three_letter_code]
        return f"AminoAcid: {self.__one_letter_code} {self.__three_letter_code} {self.is_assigned()}\nAtoms: {self.__atoms_assignment_states}\nChemical Shifts: {df_atoms}\n"

    @property
    def one_letter_code(self) -> str:
        return self.__one_letter_code

    @property
    def three_letter_code(self) -> str:
        return self.__three_letter_code


    def get_atoms(self) -> list:
        return list(self.__atoms_assignment_states.keys())


    def assign_atom(self, atom) -> None:
        match self.__atoms_assignment_states.get(atom):
            case False: self.__atoms_assignment_states[atom] = True
            case True: return
            case _: raise ValueError(f"The atom of {atom} do NOT exist! ")


    def is_assigned(self, atom = None) -> bool:
        return all(self.__atoms_assignment_states.values()) if atom is None else self.__atoms_assignment_states[atom]


class Protein:
    def __init__(self, sequence: str) -> None:
        self.__sequence = sequence
        self.__amino_acids = [AminoAcid(id) for id in sequence]


    def __str__(self) -> str:
        str_amino_acids = [f"\n{index} {aa}" for index, aa in enumerate(self.__amino_acids)]
        return f"Protein Sequence: Len = {len(self.__sequence)} {self.__sequence}\n{''.join(str_amino_acids)}"


    def get_sequence(self) -> str:
        return self.__sequence


    def __getitem__(self, index):
        return self.__amino_acids[index]


    def is_assigned(self, index = None, atom = None) -> bool:
        return all( [aa.is_assigned(atom) for aa in self.__amino_acids] ) if index is None else self.__amino_acids[index].is_assigned(atom)


    def assign(self, index, atom) -> None:
        self.__amino_acids[index].assign_atom(atom)


class Peptide(Protein):
    def __init__(self, sequence: str, range: tuple) -> None:
        super().__init__(sequence)
        self.__range = range


    def __str__(self) -> str:
        sequence = super().get_sequence()
        sequence_peptide = super().get_sequence()[self.__range[0]:self.__range[1]]
        str_amino_acids = [f"{index} {self[index]}" for index in range(self.__range[0], self.__range[1])]
        return f"Protein Sequence: Len = {len(sequence)} {sequence}\n" +\
            f"Peptide Sequence: Len = {len(sequence_peptide)} {self.__range} {sequence_peptide}\n\n" +\
            '\n'.join(str_amino_acids)


    def get_sequence(self, full_sequence = False) -> str:
        sequence = super().get_sequence()
        return sequence if full_sequence else sequence[self.__range[0]:self.__range[1]]


    def is_assigned(self, index = None, atom = None) -> bool:
        if index is not None:
            return super().is_assigned(index, atom)
        return all( [self[index].is_assigned(atom) for index in range(self.__range[0], self.__range[1])] )


    def assign(self, index, atom) -> None:
        if index < self.__range[0] or index >= self.__range[1]:
            raise IndexError( f"The parameter of index ({index}) is out of range {self.__range}!" )
        super().assign(index, atom)
