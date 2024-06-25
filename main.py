from protein import *
from assignment import *

##########################################################################################

if __name__ == "__main__":
    tau_seq = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"

    tau = Protein(tau_seq)
    core = Peptide(tau_seq, (263, 399))

    CC = PeakList('CC', ('13C', '13C'))
    NCA = PeakList('NCA', ('13C', '15N'))
    NCACB = PeakList('NCACB', ('13C', '15N'))
    NCACX = PeakList('NCACX', ('13C', '13C', '15N'))
    NCOCX = PeakList('NCOCX', ('13C', '13C', '15N'))
    CONCA = PeakList('CONCA', ('13C', '15N', '13C'))


    # print(NCACB)
    # print(NCACB.get_assigned_peaks())
    # print(NCACB.get_unassigned_peaks())

    assign(NCACB, (60.08067,127.90678), 0, core, 319, 'Thr', 'CA')
    assign(NCACB, (60.08067,127.90678), 1, core, 319, 'Thr', 'N')
    assign(NCACB, (67.47906,127.83724), 0, core, 319, 'Thr', 'CB')
    assign(NCACB, (67.47906,127.83724), 1, core, 319, 'Thr', 'N')

    assign(CONCA, (59.08024,126.09000,172.96581), 2, core, 318, 'Val', 'C')
    # assign(NCACB, (46.75326, 115.61445), 0, core, 300, 'Ala', 'CA')

    print(core[319])
    extend_assignment_left(core, 319)
