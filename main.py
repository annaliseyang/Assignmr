import pandas as pd

tau = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
core = tau[263:399] # sequence of the tau rigid core
print("Core sequence:", core)

with open('chemical_shifts.csv', 'r') as f:
    cs = pd.read_csv(f)
    print(cs.head(5))

# dictionary mapping amino acid one-letter codes to three-letter codes
aa_dict = {
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
    'V': 'Val'
}


# for i, row in cs.iterrows():
#     id = row['comp_id']
#     new_id = id[0] + str.lower(id[1:])
#     print(new_id)
#     cs['comp_id'][i] = new_id

cs.to_csv('chemical_shifts.csv', index=False)


class AminoAcid:
    def __init__(self, one_letter_code: str, atoms: list, chemical_shifts: pd.DataFrame):
        self.one_letter_code = one_letter_code
        self.three_letter_code = aa_dict[one_letter_code]
        self.atoms = atoms
        self.atoms_assignment_states = {
            atom: False for atom in atoms
        }
        self.chemical_shifts = chemical_shifts
        self.is_assigned = False

    def __str__(self):
        return self.one_letter_code

    def assign_atom(self, atom):
        self.atoms_assignment_states[atom] = True



class Protein:
    def __init__(self, sequence):
        self.sequence = sequence

    def __str__(self):
        return self.sequence

    def __get__(self, index):
        return self.sequence[index]
