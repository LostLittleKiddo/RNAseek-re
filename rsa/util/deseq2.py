# rsa/util/deseq2.py
import pandas as pd
import logging
import re
import os
from django.conf import settings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def parse_gtf_for_symbols(gtf_path):
    """Parse GTF file to extract gene_id to gene_name mapping."""
    try:
        # Read GTF, focusing on 'gene' features
        gtf_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        gtf_df = pd.read_csv(gtf_path, sep='\t', comment='#', names=gtf_cols, low_memory=False)
        gtf_df = gtf_df[gtf_df['feature'] == 'gene']
        
        # Extract gene_id and gene_name from attributes
        gene_mapping = {}
        for attr in gtf_df['attributes']:
            gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
            if gene_id_match and gene_name_match:
                gene_mapping[gene_id_match.group(1)] = gene_name_match.group(1)
        
        logger.info(f"Parsed {len(gene_mapping)} gene symbols from GTF")
        return gene_mapping
    except Exception as e:
        logger.error(f"Failed to parse GTF file {gtf_path}: {str(e)}")
        raise RuntimeError(f"Failed to parse GTF file: {str(e)}")

def prepare_metadata(meta_data, counts_data):
    """Prepare metadata and align it with counts data (from deseq2.py)."""
    meta_index = meta_data.columns[0]
    meta_data = meta_data.set_index(meta_index)
    shared_samples = counts_data.index.intersection(meta_data.index)
    if len(shared_samples) == 0:
        raise ValueError("No common samples between counts and metadata. Check sample names.")
    logger.info(f"Shared samples: {len(shared_samples)} out of {len(counts_data.index)}")
    counts_data = counts_data.loc[shared_samples]
    meta_data = meta_data.loc[shared_samples]
    return meta_data, counts_data

def run_deseq2(project, counts_file, metadata_file, output_dir):
    """
    Run DESeq2 on counts.csv and metadata.csv, adding gene symbols from GTF.

    Args:
        project: Project instance (contains species).
        counts_file: Path to counts.csv from FeatureCounts.
        metadata_file: Path to metadata.csv.
        output_dir: Directory for DESeq2 output (deseq2_results.csv).

    Returns:
        list: Path to deseq2_results.csv (single file as a list for consistency).
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "deseq2_results.csv")
    
    # Map species to GTF (reuse from featurecounts.py)
    species_to_gtf = {
        'human': 'Homo_sapiens.GRCh38.gtf',
        'mouse': 'Mus_musculus.GRCm39.gtf',
        'yeast': 'Saccharomyces_cerevisiae.R64-1-1.gtf'
    }
    gtf_file = species_to_gtf.get(project.species.lower(), None)
    if not gtf_file:
        raise RuntimeError(f"No GTF file defined for species: {project.species}")
    gtf_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'featurecounts', gtf_file)
    if not os.path.exists(gtf_path):
        raise RuntimeError(f"GTF file not found: {gtf_path}")

    # Read counts and metadata
    try:
        counts = pd.read_csv(counts_file, sep='\t')
        counts = counts.set_index('Geneid')
        columns_to_remove = ['Chr', 'Start', 'End', 'Strand', 'Length']
        counts = counts.drop(columns=columns_to_remove, errors='ignore')
        numeric_counts = counts.select_dtypes(include=['int64', 'float64'])
        counts = counts[numeric_counts.sum(axis=1) > 0]
        counts = counts.T  # Samples as rows, genes as columns

        metadata = pd.read_csv(metadata_file)
        metadata, counts = prepare_metadata(metadata, counts)
        
        # Run DESeq2
        dds = DeseqDataSet(counts=counts, metadata=metadata, design='condition')
        dds.deseq2()
        stat_res = DeseqStats(dds, n_cpus=8, contrast=['condition', metadata['condition'].unique()[0], metadata['condition'].unique()[1]])
        stat_res.summary()
        results_df = stat_res.results_df

        # Add gene symbols
        gene_mapping = parse_gtf_for_symbols(gtf_path)
        results_df['gene_symbol'] = results_df.index.map(gene_mapping)
        results_df['gene_symbol'] = results_df['gene_symbol'].fillna('Na')
        
        # Save results
        results_df.to_csv(output_file)
        logger.info(f"DESeq2 results saved to: {output_file}")

        # Register output file
        file_size = os.path.getsize(output_file) if os.path.isfile(output_file) else None
        ProjectFiles.objects.create(
            project=project,
            type='deseq2_results',
            path=output_file,
            is_directory=False,
            file_format='csv',
            size=file_size
        )
        logger.info(f"Registered DESeq2 output: {output_file} with size {file_size} bytes")

        return [output_file]
    
    except Exception as e:
        logger.error(f"DESeq2 failed: {str(e)}")
        raise RuntimeError(f"DESeq2 failed: {str(e)}")