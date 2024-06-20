from protein import *
from assignment import *

##########################################################################################

if __name__ == "__main__":
    tau_seq = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
    # core = tau[263:399] # sequence of the tau rigid core

    # with open('chemical_shifts.csv', 'r') as f:
    #     cs = pd.read_csv(f)
    # #    print(cs.head(5))

    #     aa_1 = AminoAcid('A')
    #     aa_1.assign_atom(aa_1.get_atoms()[1])
    #     print(aa_1.one_letter_code)
    #     print(aa_1.three_letter_code)
    #     print(aa_1.is_assigned(aa_1.get_atoms()[0]))
    #     print(aa_1.is_assigned(aa_1.get_atoms()[1]))
    #     print(aa_1.is_assigned())
    #     print(aa_1.get_atoms())
    #     print(aa_1)

    #     aa_2 = AminoAcid('Glu')
    #     aa_2.assign_atom(aa_2.get_atoms()[1])
    #     aa_2.assign_atom(aa_2.get_atoms()[2])
    #     aa_2.assign_atom(aa_2.get_atoms()[3])
    #     print(aa_2)

    #     protein_1 = Protein(tau)
    #     protein_1.assign(10, protein_1[10].get_atoms()[1])
    #     print(protein_1)
    #     print(protein_1.get_sequence())
    #     print(protein_1.is_assigned(1, protein_1[1].get_atoms()[1]))
    #     print(protein_1.is_assigned(1))
    #     print(protein_1.is_assigned())
    #     print(protein_1[0])

    #     peptide_1 = Peptide(tau, (10, 15))
    #     peptide_1.assign(10, peptide_1[10].get_atoms()[1])
    #     print(peptide_1)
    #     print(peptide_1.get_sequence())
    #     print(peptide_1.is_assigned(1, peptide_1[1].get_atoms()[1]))
    #     print(peptide_1.is_assigned(1))
    #     print(peptide_1.is_assigned())
    #     print(peptide_1[0])
    tau = Protein(tau_seq)
    core = Peptide(tau_seq, (263, 399))
    NCACB = PeakList('NCACB', ('13C', '15N'))

    print(NCACB)
    print(NCACB.get_assigned_peaks())
    assign(NCACB, (46.75326, 115.61445), '13C', core, 300, 'Ala', 'CA')

    print(core[300])
