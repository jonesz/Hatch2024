import json

def check_accession_presence(rs_data, rs_ID, accession_ID):

    # Search for rs_ID in the JSON data.
    variant_info = None
    for key in rs_data['primary_snapshot_data']['ga4gh']:
        if key.split('.')[0] == accession_ID:
            variant_info = key
            break

    if variant_info:
        print(f"{rs_ID} variant IS present in {accession_ID} sequence.")
        return True
    else:
        print(f"Variant {rs_ID} was not found in the downloaded data.")
        return False

    # # Check if accesion_ID is associated with snp_ID.
    # if variant_info:
    #     if accession_ID in variant_info['associated_sequences']:
    #         print(f"{rs_ID} variant IS present in {accession_ID} sequence.")
    #         return True
    #     else:
    #         print(f"{rs_ID} variant IS present in {accession_ID} sequence.")
    #         return False
    # else:
    #     print(f"Variant {rs_ID} was not found in the downloaded data.")
    #     return False

def main():

    rs_filepath = "dbsnp/1801133.json"
    rs_ID = "1801133"

    astro_accessions_filepath = "blast/astro-accesions.json"

    results_dict = {"rs1801133":{}}
    results_filepath = "astro-snp-results.json"

    with open(rs_filepath, 'r') as file:
        rs_data = json.load(file)

    with open(astro_accessions_filepath, 'r') as file:
        astro_accession_IDs = json.load(file)

    for astro, accession_ID in astro_accession_IDs.items():
        result = check_accession_presence(rs_data, rs_ID, accession_ID)
        print(f"{astro} has rs_ID: {result}")
        results_dict["rs1801133"][astro] = result

    results_string = json.dumps(results_dict)
    with open(results_filepath, "w") as file:
        file.write(results_string)

if __name__ == "__main__":
    main()