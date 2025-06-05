import pandas as pd
import logging
from pydeseq2.dds import DeseqDataSet

# Configure basic logging
logging.basicConfig(level=logging.INFO)

def prepare_metadata(meta_data, counts_data):
    """Prepare metadata and align it with counts data."""
    meta_index = meta_data.columns[0]
    meta_data = meta_data.set_index(meta_index)
    shared_samples = counts_data.index.intersection(meta_data.index)
    if len(shared_samples) == 0:
        raise ValueError("No common samples between counts and metadata. Check sample names.")
    logging.info(f"Shared samples: {len(shared_samples)} out of {len(counts_data.index)}")
    counts_data = counts_data.loc[shared_samples]
    meta_data = meta_data.loc[shared_samples]
    return meta_data, counts_data

# Read the counts.csv file
counts = pd.read_csv('counts.csv', sep='\t')

# Set 'Geneid' as the index
counts = counts.set_index('Geneid')

# Remove unwanted columns
columns_to_remove = ['Chr', 'Start', 'End', 'Strand', 'Length']
counts = counts.drop(columns=columns_to_remove, errors='ignore')

# Select only numeric columns for the sum operation
numeric_counts = counts.select_dtypes(include=['int64', 'float64'])

# Filter rows where the sum of numeric columns is greater than 0
counts = counts[numeric_counts.sum(axis=1) > 0]

# Transpose counts so samples are rows and genes are columns
counts = counts.T

# Read metadata.csv
metadata = pd.read_csv('metadata.csv')

# Prepare metadata and align with counts
metadata, counts = prepare_metadata(metadata, counts)

# Print the 'condition' column from metadata
print("Condition column:", metadata['condition'].to_list())

print("Metadata columns:", metadata.columns)

dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design='~1',  # Use intercept-only model to avoid replicate issue
)

dds.deseq2()

print(dds.obs)