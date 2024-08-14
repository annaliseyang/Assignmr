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
    NCOCX = PeakList('NCOCX', ('13C', '13C', '15N')) # Cx, CO, N
    CONCA = PeakList('CONCA', ('13C', '15N', '13C')) # Ca, N, CO


    # print(NCACB)
    # print(NCACB.get_assigned_peaks())
    # print(NCACB.get_unassigned_peaks())

    assign(NCACB, (60.08067,127.90678), 0, core, 319, 'Thr', 'CA')
    assign(NCACB, (60.08067,127.90678), 1, core, 319, 'Thr', 'N')
    assign(NCACB, (67.47906,127.83724), 0, core, 319, 'Thr', 'CB')
    assign(NCACB, (67.47906,127.83724), 1, core, 319, 'Thr', 'N')

    core[305].assign_atom('C', Assignment(172.812, True))
    core[305].assign_atom('N', Assignment(112.32, True))
    core[305].assign_atom('CA', Assignment(52.329, True))
    core[305].assign_atom('CB', Assignment(63.353, True))

    start_index = 319

    print(core[start_index])
    extend_assignment_left(core, start_index)
    complete_assignment(core, start_index - 1)

    # assign(CONCA, (59.08024,126.09000,172.96581), 2, core, 318, 'Val', 'C')
    # assign(NCACB, (46.75326, 115.61445), 0, core, 300, 'Ala', 'CA')

    # print(core[319])
    # index = 319
    # while index > 310:
    #     extend_assignment_left(core, index)
    #     complete_assignment(core, index - 1)
    #     index -= 1
    #     break # break for testing

    print(core[318])
    # print(core)

    # assign(NCACB, (58.82997,126.11562), 0, core, 318, 'Val', 'CA')
    # assign(NCACB, (58.82997,126.11562), 1, core, 318, 'Val', 'N')

    # print(core[318])
    # print(core[319])

    # extend_assignment_left(core, 318)
    # NCACB.plot_peaks_interactive()
    # CONCA.plot_peaks_interactive_3d()
