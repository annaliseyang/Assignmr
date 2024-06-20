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


    def get_chemical_shifts(amino_acid = None):
        if amino_acid is None:
            return AminoAcid.__chemical_shifts
        match(len(amino_acid)):
            case 1:
                amino_acid = AminoAcid.__amino_acid_dict[amino_acid] 
                return AminoAcid.__chemical_shifts[AminoAcid.__chemical_shifts['comp_id'] == amino_acid]
            case 3:
                return AminoAcid.__chemical_shifts[AminoAcid.__chemical_shifts['comp_id'] == amino_acid]
        return ValueError( f"The parameter of amino_acid ({amino_acid}) is invalid!" )


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
        atoms_str = [ f"\033[1;33m{k}\033[0m" if v else f"{k}" for k, v in self.__atoms_assignment_states.items() ]
        return f"AminoAcid: {self.__one_letter_code} {self.__three_letter_code} {self.is_assigned()}\nAtoms: ({len(atoms_str)}) {", ".join(atoms_str)}\n"

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


if __name__ == "__main__":
    tau = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
    core = tau[263:399] # sequence of the tau rigid core

    with open('chemical_shifts.csv', 'r') as f:
        cs = pd.read_csv(f)
    #    print(cs.head(5))

        aa_1 = AminoAcid('A')
        aa_1.assign_atom(aa_1.get_atoms()[1])
        print(aa_1.one_letter_code)
        print(aa_1.three_letter_code)
        print(aa_1.is_assigned(aa_1.get_atoms()[0]))
        print(aa_1.is_assigned(aa_1.get_atoms()[1]))
        print(aa_1.is_assigned())
        print(aa_1.get_atoms())
        print(aa_1)

        aa_2 = AminoAcid('Glu')
        aa_2.assign_atom(aa_2.get_atoms()[1])
        aa_2.assign_atom(aa_2.get_atoms()[2])
        aa_2.assign_atom(aa_2.get_atoms()[3])
        print(aa_2)

        protein_1 = Protein(tau)
        protein_1.assign(10, protein_1[10].get_atoms()[1])
        # print(protein_1)
        print(protein_1.get_sequence())
        print(protein_1.is_assigned(1, protein_1[1].get_atoms()[1]))
        print(protein_1.is_assigned(1))
        print(protein_1.is_assigned())
        print(protein_1[0])

        peptide_1 = Peptide(tau, (10, 15))
        peptide_1.assign(10, peptide_1[10].get_atoms()[1])
        print(peptide_1)
        print(peptide_1.get_sequence())
        print(peptide_1.is_assigned(1, peptide_1[1].get_atoms()[1]))
        print(peptide_1.is_assigned(1))
        print(peptide_1.is_assigned())
        print(peptide_1[0])

        print("=======================================")
        print(AminoAcid.get_chemical_shifts())
        print(AminoAcid.get_chemical_shifts('P'))
        print(AminoAcid.get_chemical_shifts('Pro'))
