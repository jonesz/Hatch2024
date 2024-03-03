from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Load your FASTA file
astro_MTHFRs = (
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro1_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro2_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro3_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro4_MTHFR.fasta",
    "astro-nutrition-dataset/MTHFR Astro_Revised/Astro5_MTHFR.fasta"
)

# Read the FASTA file
# sequences = SeqIO.parse(astro_MTHFRs[0], "fasta")
sequences = SeqIO.parse(astro_MTHFRs[0], "fasta")

# Perform BLAST search for each sequence
for seq_record in sequences:
    result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.seq)
    print(result_handle)

    # Check status of the request
    # Note: This will print "Status: IN PROGRESS" until the request is completed

    # Save the BLAST results to a file
    output_file = f"{seq_record.id}_blast_results.xml"
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()

print("Done!")