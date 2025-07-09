import csv
import re
import os

def is_kegg_id(term):
    """Check if the term matches KEGG ID format (hsaXXXXX, NXXXXX, athXXXXX, osaXXXXX, or zmaXXXXX)."""
    return bool(re.match(r'^(hsa|N|ath|osa|zma)\d{5}$', term))

def is_go_or_hp_id(term):
    """Check if the term matches GO:XXXXXXX or HP:XXXXXXX format (7 digits)."""
    return bool(re.match(r'^(GO|HP):\d{7}$', term))

def process_gmt_file(input_file, output_file):
    """Process GMT file and create CSV with term name and corresponding link."""
    with open(input_file, 'r') as gmt_file, open(output_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Term_Name', 'Link'])  # Write header

        for line in gmt_file:
            columns = line.strip().split('\t')
            if len(columns) < 2:
                continue  # Skip invalid lines

            # Check both possible columns for KEGG or GO/HP ID
            term_id = None
            term_name = None
            if is_kegg_id(columns[0]) or is_go_or_hp_id(columns[0]):
                term_id = columns[0]
                term_name = columns[1]
            elif is_kegg_id(columns[1]) or is_go_or_hp_id(columns[1]):
                term_id = columns[1]
                term_name = columns[0]
            else:
                continue  # Skip if no valid ID found

            # Create appropriate link based on ID prefix
            if term_id.startswith('hsa'):
                link = f"https://www.genome.jp/dbget-bin/www_bget?pathway:{term_id}"
            elif term_id.startswith('N'):
                link = f"https://www.genome.jp/entry/{term_id}"
            elif term_id.startswith('ath'):
                link = f"https://www.genome.jp/dbget-bin/www_bget?pathway:{term_id}"
            elif term_id.startswith('osa'):
                link = f"https://www.genome.jp/dbget-bin/www_bget?pathway:{term_id}"
            elif term_id.startswith('zma'):
                link = f"https://www.genome.jp/dbget-bin/www_bget?pathway:{term_id}"
            elif term_id.startswith('GO:'):
                link = f"https://www.ebi.ac.uk/QuickGO/term/{term_id}"
            elif term_id.startswith('HP:'):
                link = f"https://next.monarchinitiative.org/{term_id}"
            else:
                continue  # Skip if ID is not recognized

            csv_writer.writerow([term_name, link])

if __name__ == "__main__":
    input_directory = os.path.dirname(__file__)
    output_directory = os.path.join(input_directory, "converted_links")

    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Process each GMT file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".gmt"):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}_links.csv")
            process_gmt_file(input_file, output_file)