from Bio import Entrez, SeqIO
import json
import os

email = os.environ.get('EMAIL')
print(email)

# Define your email address for Entrez
Entrez.email = os.environ.get('EMAIL')

# Function to fetch sequence by accession number
def fetch_sequence(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record.seq

def check_snp_presence(accession, snp_id, ref_nucleotide, snp_position):
    # Fetch the sequence using Entrez
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    sequence = str(record.seq)

    # Check if the SNP position matches the reference sequence
    if sequence[snp_position - 1] == ref_nucleotide:
        print(f"The SNP {snp_id} is present at position {snp_position}.")
        return True
    else:
        print(f"The SNP {snp_id} is not present at position {snp_position}.")
        return False

# Main code
def main():

    acc_filepath = "BLAST/astro-accesions.json"

    with open(acc_filepath, 'r') as file:
        acc_dict = json.load(file)

    astro_IDs = sorted(acc_dict.keys())

    # rsID for the SNP of interest
    snp_id = "rs1801133"
    ref_nucleotide = 'C'
    snp_position = 677

    results_dict = {"rs1801133": {}}
    results_filepath = "astro-snp-results.json"

    print("\n\n")

    for astro_ID in astro_IDs:

        # Fetch sequence
        sequence = fetch_sequence(acc_dict[astro_ID])

        result = check_snp_presence(acc_dict[astro_ID], snp_id, ref_nucleotide, snp_position)
        print(f"{astro_ID} has rs1801133 variant: {result}\n")

        results_dict["rs1801133"][astro_ID] = result

    results_string = json.dumps(results_dict)
    
    with open(results_filepath, "w") as file:
        file.write(results_string)

if __name__ == "__main__":

    main()
