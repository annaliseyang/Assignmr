from protein import *
from assignment import *
import time

##########################################################################################

def main(protein, start_index, save):
    """
    Main function to recursively assign chemical shifts to the protein.
    start index: a fully assigned amino acid where recursive assignments will start from.
    """
    start_time = time.time()
    try:
        print(protein[start_index])
        extend_assignment_left(protein, start_index)
        complete_assignment(protein, start_index - 1)
    except Exception as e:
        print(f"Error extending assignment: {e}")
        protein.to_pickle(save + '.pkl')
        protein.to_txt(save + '.txt')
    finally:
        run_time = time.time() - start_time
        print(f"Run time: {run_time//60} minutes {run_time%60} seconds")




if __name__ == "__main__":

    tau_seq = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"

    tau = Protein(tau_seq)
    i_min, i_max = (263, 399)
    core = Peptide(tau_seq, (i_min, i_max))

    CC = PeakList('CC', ('13C', '13C'))
    NCA = PeakList('NCA', ('13C', '15N'))
    NCACB = PeakList('NCACB', ('13C', '15N'))
    NCACX = PeakList('NCACX', ('13C', '13C', '15N'))
    NCOCX = PeakList('NCOCX', ('13C', '13C', '15N')) # Cx, CO, N
    CONCA = PeakList('CONCA', ('13C', '15N', '13C')) # Ca, N, CO


    # print(NCACB)
    # print(NCACB.get_assigned_peaks())
    # print(NCACB.get_unassigned_peaks())
    thr_indices = [i for i in core.get_indices_of_amino_acid('Thr') if i in range(i_min, i_max+1)]
    print(f"Thr indices: {thr_indices}") # [263, 319, 361, 373, 377, 386]

    thr1_index = 386
    core[thr1_index].assign_atom('C', Assignment(171.151, True))
    core[thr1_index].assign_atom('N', Assignment(127.912, True))
    core[thr1_index].assign_atom('CA', Assignment(60.063, True))
    core[thr1_index].assign_atom('CB', Assignment(67.419, True))
    # assign(NCACB, (60.08067,127.90678), 0, core, 319, 'Thr', 'CA')
    # assign(NCACB, (60.08067,127.90678), 1, core, 319, 'Thr', 'N')
    # assign(NCACB, (67.47906,127.83724), 0, core, 319, 'Thr', 'CB')
    # assign(NCACB, (67.47906,127.83724), 1, core, 319, 'Thr', 'N')

    core[305].assign_atom('C', Assignment(172.812, True))
    core[305].assign_atom('N', Assignment(112.32, True))
    core[305].assign_atom('CA', Assignment(52.329, True))
    core[305].assign_atom('CB', Assignment(63.353, True))

    lys_indices = [i for i in core.get_indices_of_amino_acid('Lys') if i in range(i_min, i_max+1)]
    print(f"Lys indices: {lys_indices}") # [267, 274, 280, 281, 290, 294, 298, 311, 317, 321, 331, 340, 343, 347, 353, 369, 370, 375, 383, 385, 395]
    lys1_index = 290
    core[lys1_index].assign_atom('CA', Assignment(56.365, True))
    core[lys1_index].assign_atom('N', Assignment(127.471, True))
    core[lys1_index].assign_atom('CB', Assignment(29.894, True))

    start_index = thr1_index

    main(core, start_index, 'core')


    # assign(CONCA, (59.08024,126.09000,172.96581), 2, core, 318, 'Val', 'C')
    # assign(NCACB, (46.75326, 115.61445), 0, core, 300, 'Ala', 'CA')

    # print(core[319])
    # index = 319
    # while index > 310:
    #     extend_assignment_left(core, index)
    #     complete_assignment(core, index - 1)
    #     index -= 1
    #     break # break for testing

    # print(core[318])
    # print(core)

    # assign(NCACB, (58.82997,126.11562), 0, core, 318, 'Val', 'CA')
    # assign(NCACB, (58.82997,126.11562), 1, core, 318, 'Val', 'N')

    # print(core[318])
    # print(core[319])

    # extend_assignment_left(core, 318)
    # NCACB.plot_peaks_interactive()
    # CONCA.plot_peaks_interactive_3d()
