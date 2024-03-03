from Bio import pairwise2, SeqIO
from Bio.Align import PairwiseAligner
from Bio.pairwise2 import format_alignment

astro_MTHFR = (
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro1_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro2_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro3_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro4_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro5_MTHFR.fasta"
)
ref_MTHFR = "astro-nutrition-dataset/MTHFR Gene Comp/MTHFR_Gene.fasta"

chr_1 = "astro-nutrition-dataset/HomoSapien_Chr1.fasta"

# Read the FASTA sequences from files
seq1_record = SeqIO.read(astro_MTHFR[0], "fasta")
# seq1_record = SeqIO.read(ref_MTHFR, "fasta")
seq2_record = SeqIO.read(astro_MTHFR[2], "fasta")

# Convert sequences to Biopython Seq objects
seq1 = str(seq1_record.seq)
seq2 = str(seq2_record.seq)

alignments = pairwise2.align.globalxx(seq1, seq2)
best_alignment = max(alignments, key=lambda x: x[2])

# Write alignment results to a new file
with open("alignment_results.txt", "w") as f:
    f.write(format_alignment(*best_alignment, full_sequences=True))

print("Done! Alignment written to file.")