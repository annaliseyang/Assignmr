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


    def get_chemical_shifts(amino_acid = None, atom = None):
        if amino_acid is None:
            return AminoAcid.__chemical_shifts
        if atom:
            return AminoAcid.__chemical_shifts[(AminoAcid.__chemical_shifts['comp_id'] == amino_acid) & (AminoAcid.__chemical_shifts['atom_id'] == atom)]
        match(len(amino_acid)):
            case 1:
                amino_acid = AminoAcid.__amino_acid_dict[amino_acid]
                return AminoAcid.__chemical_shifts[AminoAcid.__chemical_shifts['comp_id'] == amino_acid]
            case 3:
                return AminoAcid.__chemical_shifts[AminoAcid.__chemical_shifts['comp_id'] == amino_acid]
        return ValueError( f"The parameter of amino_acid ({amino_acid}) is invalid!" )


    def get_chemical_shift_range(amino_acid, atom, num_std = 2) -> tuple:
        cs = AminoAcid.get_chemical_shifts(amino_acid, atom)
        avg = float(cs['avg'].iloc[0])
        half_range = float(cs['std'].iloc[0]) * num_std
        return (avg - half_range, avg + half_range)

    def get_cx_chemical_shift_ranges(amino_acid, num_std = 2) -> dict:
        '''
        Return a dictionary of the chemical shift ranges of all the Cx atoms of the amino acid, excluding C.
        '''
        ranges = {
            atom: AminoAcid.get_chemical_shift_range(amino_acid, atom, num_std)
            for atom in AminoAcid.get_chemical_shifts(amino_acid)['atom_id']
            if atom.startswith('C') and atom != 'C'
        }
        # for atom, shift in AminoAcid.get_chemical_shifts(amino_acid).iterrows():
        #     if atom.startswith('C') and atom != 'C':
        #         print(atom, shift)
        # ca_range = AminoAcid.get_chemical_shift_range(amino_acid, 'CA', num_std)
        # cb_range = AminoAcid.get_chemical_shift_range(amino_acid, 'CB', num_std)
        return ranges

    def __init__(self, id: str):
        id_2 = self.__amino_acid_dict[id]
        match(len(id)):
            case 1:
                self.__one_letter_code = id
                self.__three_letter_code = id_2
            case 3:
                self.__one_letter_code = id_2
                self.__three_letter_code = id
        cs = AminoAcid.__chemical_shifts

        condition = [
                (cs['atom_id'][row].startswith('C') or cs['atom_id'][row].startswith('N'))
                and cs['comp_id'][row] == self.__three_letter_code
                for row in range(len(cs['atom_id']))
            ]

        df_atoms = cs[condition]
        self.__atoms_assignments = { atom: False for atom in df_atoms['atom_id'].to_list() }


    def __str__(self) -> str:
        atoms_str = [ f"\033[1;33m{k}\033[0m" if v else f"{k}" for k, v in self.__atoms_assignments.items() ]
        return f"{self.__three_letter_code} ({self.one_letter_code}) {self.is_assigned()}\nAtoms: ({len(atoms_str)}) {', '.join(atoms_str)}\n"

    @property
    def one_letter_code(self) -> str:
        return self.__one_letter_code

    @property
    def three_letter_code(self) -> str:
        return self.__three_letter_code


    def get_atoms(self) -> list:
        return list(self.__atoms_assignments.keys())


    def assign_atom(self, atom: str, assignment) -> None:
        if self.__atoms_assignments.get(atom) == None:
            raise ValueError(f"The atom of {atom} do NOT exist! ")
        self.__atoms_assignments[atom] = assignment


    def get_assignment(self, atom: str):
        assignment = self.__atoms_assignments[atom]
        return assignment if assignment != False else None


    def is_assigned(self, atom: str = None) -> bool:
        if atom is None:
            return all(assignment != False for assignment in self.__atoms_assignments.values())
        return self.__atoms_assignments[atom] != False


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
        return self.__amino_acids[index-1] # input index starts from 1


    def assign_atom(self, index: int, atom: str , assignment) -> None:
        self[index].assign_atom(atom, assignment)


    def is_assigned(self, index: int = None, atom: str = None) -> bool:
        return all( [aa.is_assigned(atom) for aa in self.__amino_acids] ) if index is None else self.__amino_acids[index].is_assigned(atom)


    def get_assignment(self, index: int, atom: str):
        return self[index].get_assignment(atom)


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


    def assign_atom(self, index: int, atom: str , assignment) -> None:
        if index < self.__range[0] or index >= self.__range[1]:
            raise IndexError( f"The parameter of index ({index}) is out of range {self.__range}!" )
        super().assign_atom(index, atom, assignment)


if __name__ == "__main__":
    tau = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
    core = tau[263:399] # sequence of the tau rigid core

    assignment_1 = (123.45, True)
    assignment_2 = (789.00, False)

    # print(assignment_1, assignment_2)
    # print(assignment_1.chemical_shift, assignment_2.tentative)


    with open('chemical_shifts.csv', 'r') as f:
        cs = pd.read_csv(f)
    #    print(cs.head(5))

        aa_1 = AminoAcid('A')
        # aa_1.assign_atom(aa_1.get_atoms()[0])
        # print(aa_1.one_letter_code)
        # print(aa_1.three_letter_code)
        # print(aa_1.is_assigned(aa_1.get_atoms()[0]))
        # print(aa_1.is_assigned(aa_1.get_atoms()[1]))
        # print(aa_1.is_assigned())
        # print(aa_1.get_atoms())
        print(aa_1)
        print("Assignment: ", aa_1.get_assignment(aa_1.get_atoms()[0]))

        aa_2 = AminoAcid('Glu')
        aa_2.assign_atom(aa_2.get_atoms()[1], assignment_1)
        aa_2.assign_atom(aa_2.get_atoms()[2], assignment_1)
        aa_2.assign_atom(aa_2.get_atoms()[3], assignment_1)
        print(aa_2)
        print("Assignment: ", aa_2.get_assignment(aa_2.get_atoms()[0]))
        print("Assignment: ", aa_2.get_assignment(aa_2.get_atoms()[3]))

        protein_1 = Protein(tau)
        protein_1.assign_atom(10, protein_1[10].get_atoms()[1], assignment_2)
        # print(protein_1)
        print(protein_1.get_sequence())
        print(protein_1.is_assigned(1, protein_1[1].get_atoms()[1]))
        print(protein_1.is_assigned(1))
        print(protein_1.is_assigned())
        print(protein_1[0])

        peptide_1 = Peptide(tau, (10, 15))
        peptide_1.assign_atom(10, peptide_1[10].get_atoms()[1], assignment_2)
        print(peptide_1)
        print(peptide_1.get_sequence())
        print(peptide_1.is_assigned(1, peptide_1[1].get_atoms()[1]))
        print(peptide_1.is_assigned(1))
        print(peptide_1.is_assigned())
        print(peptide_1[0])
        print("Assignment: ", peptide_1.get_assignment(10, peptide_1[10].get_atoms()[1]))
        print("Assignment: ", peptide_1.get_assignment(11, peptide_1[11].get_atoms()[1]))

        print("=======================================")
        print(AminoAcid.get_chemical_shifts())
        print(AminoAcid.get_chemical_shifts('P'))
        print(AminoAcid.get_chemical_shifts('Pro'))
