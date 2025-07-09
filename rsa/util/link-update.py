import pandas as pd

# Read the input CSV files

species_to_link = {
    'human': {
        'go': 'homo_sapiens_go_links.csv',
        'kegg': 'homo_sapiens_kegg_links.csv'
    },
    'mouse': {
        'go': 'mus_musculus_go_links.csv',
        'kegg': 'mus_musculus_kegg_links.csv'
    },
    'yeast': {
        'go': 'saccharomyces_cerevisiae_go_links.csv',
        'kegg': 'saccharomyces_cerevisiae_kegg_links.csv'
    },
    'arabidopsis': {
        'go': 'arabidopsis_thaliana_go_links.csv',
        'kegg': 'arabidopsis_thaliana_kegg_links.csv'
    },
    'worm': {
        'go': 'caenorhabditis_elegans_go_links.csv',
        'kegg': 'caenorhabditis_elegans_kegg_links.csv'
    },
    'zebrafish': {
        'go': 'danio_rerio_go_links.csv',
        'kegg': 'danio_rerio_kegg_links.csv'
    },
    'fly': {
        'go': 'drosophila_melanogaster_go_links.csv',
        'kegg': 'drosophila_melanogaster_kegg_links.csv'
    },
    'rice': {
        'go': 'oryza_sativa_go_links.csv',
        'kegg': 'oryza_sativa_kegg_links.csv'
    },
    'maize': {
        'go': 'zea_mays_go_links.csv',
        'kegg': 'zea_mays_kegg_links.csv'
    }
}
csv1 = pd.read_csv('go_gsea_results.csv')
csv2 = pd.read_csv('saccharomyces_cerevisiae_go_links.csv')

# Merge the dataframes based on Term and Term_Name
merged_df = pd.merge(csv1, csv2[['Term_Name', 'Link']], 
                    left_on='Term', 
                    right_on='Term_Name', 
                    how='left')

# Drop the redundant Term_Name column
merged_df = merged_df.drop('Term_Name', axis=1)

# Reorder columns to place Link as the second column
columns = merged_df.columns.tolist()
term_index = columns.index('Term')
columns.insert(term_index + 1, columns.pop(columns.index('Link')))
merged_df = merged_df[columns]

# Save the result to a new CSV file
merged_df.to_csv('go_gsea_results.csv', index=False)

print("CSV files merged successfully. Output saved to 'go_gsea_results.csv'")