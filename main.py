import pandas as pd

######################################################################

class AminoAcid:
    amino_acid_dict = {
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

        'Ala': 'A',
        'Arg': 'R',
        'Asn': 'N',
        'Asp': 'D',
        'Cys': 'C',
        'Gln': 'Q',
        'Glu': 'E',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Leu': 'L',
        'Lys': 'K',
        'Met': 'M',
        'Phe': 'F',
        'Pro': 'P',
        'Ser': 'S',
        'Thr': 'T',
        'Trp': 'W',
        'Tyr': 'Y',
        'Val': 'V'
    }


    def __init__(self, letter_code: str, atoms: list, chemical_shifts: pd.DataFrame):
        letter_code_2 = self.amino_acid_dict[letter_code]

        match(len(letter_code)):
            case 1:
                self.one_letter_code = letter_code
                self.three_letter_code = letter_code_2 
            case 3:
                self.one_letter_code = letter_code_2
                self.three_letter_code = letter_code 

        self.atoms = atoms
        self.atoms_assignment_states = { atom: False for atom in atoms }
        self.chemical_shifts = chemical_shifts


    def __str__(self):
        return f"AminoAcid: {self.one_letter_code} {self.three_letter_code} {self.is_assigned_all()}\nAtoms: {self.atoms}\nChemical Shifts: {self.chemical_shifts}"


    def assign_atom(self, atom):
        self.atoms_assignment_states[atom] = True


    def get_atom_assign(self, atom):
        return self.atoms_assignment_states[atom]


    def is_assigned_all(self):
        return all(self.atoms_assignment_states.values())


class Protein:
    def __init__(self, sequence, atoms: list, chemical_shifts: pd.DataFrame):
        self.sequence = sequence
        self.atoms = atoms
        self.chemical_shifts = chemical_shifts
        self.amino_acids = [ AminoAcid( letter_code, atoms, chemical_shifts) for letter_code in sequence]

    def __str__(self):
        return f"sequence Len = {len(self.sequence)} {self.sequence}\nAtoms: {self.atoms}\nChemical Shifts: {self.chemical_shifts}"

    def __getitem__(self, index):
        return self.amino_acids[index]

##########################################################################################

tau = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
core = tau[263:399] # sequence of the tau rigid core
print("Core Sequence:", core)

with open('chemical_shifts.csv', 'r') as f:
    cs = pd.read_csv(f)
    print(cs.head(5))

    aa_1 = AminoAcid('A', ['A', 'B', 'C'], cs)
    aa_1.assign_atom('B')
    print(aa_1)

    aa_2 = AminoAcid('Glu', ['a', 'b', 'c'], cs)
    aa_2.assign_atom('a')
    aa_2.assign_atom('b')
    aa_2.assign_atom('c')
    print(aa_2)

    protein_1 = Protein(core, 'XYZ', cs)
    print(protein_1)
    print(protein_1[10])
