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


    print(NCACB)
    print(NCACB.get_assigned_peaks())
    print(NCACB.get_unassigned_peaks())
    assign(NCACB, (46.75326, 115.61445), '13C', core, 300, 'Ala', 'CA')

    print(core[300])

    # find peaks in CONCA peaklist given N and CA chemical shifts
    # peaks = CONCA.get_peaks_along_dimension((46.75326, 115.61445), 2, tolerance=10)
    Val_C_chemical_shifts = AminoAcid.get_chemical_shifts('Val', 'C')
    print(Val_C_chemical_shifts)
    target_range = AminoAcid.get_chemical_shift_range('Val', 'C')
    print('target range:', target_range)
    peaks = CONCA.get_peaks_along_dimension((60.08067,127.90678), 2, tolerance=1, num_peaks=1, cs_range=target_range)

    print_filtered(peaks)

    # NCACB.plot_peaks()
    # NCACB.plot_peaks(assigned_only=True)
    NCACB.plot_peaks_interactive()
