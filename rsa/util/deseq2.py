# rsa/util/deseq2.py
import pandas as pd
import logging
import re
import os
from django.conf import settings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from rsa.models import Project, ProjectFiles
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

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

def create_cluster_heatmap(counts, output_path, project):
    """Generate a clustered heatmap with dendrograms, styled like the provided example."""
    try:
        # Subset to top 50 most variable genes for MVP (same as before)
        counts_var = counts.var(axis=0).sort_values(ascending=False).head(50).index
        subset_counts = counts[counts_var]
        
        # Normalize data to range [-2, 2] to match the example
        min_val = subset_counts.min().min()
        max_val = subset_counts.max().max()
        range_val = max(abs(min_val), abs(max_val))
        normalized_counts = 2 * (subset_counts - min_val) / (max_val - min_val) - 1 if max_val != min_val else subset_counts
        
        # Generate heatmap with clustering
        plt.figure(figsize=(10, 8))
        sns.heatmap(normalized_counts, cmap='RdBu_r', vmin=-2, vmax=2, center=0,
                    xticklabels=True, yticklabels=True, 
                    cbar_kws={'label': 'Normalized Value'})
        plt.title(f"Clustered Heatmap - {project.name}")
        plt.xlabel("Samples")
        plt.ylabel("Genes")
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()
        logger.info(f"Clustered heatmap saved to: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Failed to create clustered heatmap: {str(e)}")
        raise RuntimeError(f"Clustered heatmap failed: {str(e)}")

def create_pca_plot(counts, metadata, output_path, project):
    """Generate a basic PCA plot of samples based on counts."""
    try:
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(counts)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Condition': metadata['condition']
        })
        
        plt.figure(figsize=(8, 6))
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Condition', s=100)
        plt.title(f"PCA Plot - {project.name}")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()
        logger.info(f"PCA plot saved to: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Failed to create PCA plot: {str(e)}")
        raise RuntimeError(f"PCA plot failed: {str(e)}")

def run_deseq2(project, counts_file, metadata_file, output_dir):
    """
    Run DESeq2 on counts.csv and metadata.csv, adding gene symbols from GTF.
    Filter results by project.pvalue_cutoff, log2FoldChange > 1, and baseMean > 10.
    Generate cluster heatmap, PCA.

    Args:
        project: Project instance (contains species and pvalue_cutoff).
        counts_file: Path to counts.csv from FeatureCounts.
        metadata_file: Path to metadata.csv.
        output_dir: Directory for DESeq2 output (deseq2_results.csv and plots).

    Returns:
        list: Paths to deseq2_results.csv and visualization PNGs.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "deseq2_results.csv")
    heatmap_output = os.path.join(output_dir, "cluster_heatmap.png")
    pca_output = os.path.join(output_dir, "pca_plot.png")

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
        stat_res = DeseqStats(dds, n_cpus=1, contrast=['condition', metadata['condition'].unique()[0], metadata['condition'].unique()[1]])
        stat_res.summary()
        results_df = stat_res.results_df

        # Add gene symbols
        gene_mapping = parse_gtf_for_symbols(gtf_path)
        results_df['gene_symbol'] = results_df.index.map(gene_mapping)
        results_df['gene_symbol'] = results_df['gene_symbol'].fillna('Na')
        
        # Filter results based on pvalue_cutoff, log2FoldChange, and baseMean
        pvalue_cutoff = project.pvalue_cutoff
        results_df = results_df[
            (results_df['padj'] < pvalue_cutoff) &
            (results_df['log2FoldChange'].abs() > 1.0) &
            (results_df['baseMean'] > 10.0)
        ]
        logger.info(f"Filtered DESeq2 results to {len(results_df)} genes with padj < {pvalue_cutoff}, "
                    f"|log2FoldChange| > 1.0, and baseMean > 10.0")
        
        # Save results
        results_df.to_csv(output_file)
        logger.info(f"DESeq2 results saved to: {output_file}")

        # Generate visualizations
        output_files = [output_file]
        create_cluster_heatmap(counts, heatmap_output, project)
        create_pca_plot(counts, metadata, pca_output, project)

        # Register visualization files
        for plot_path in [heatmap_output, pca_output]:
            if os.path.exists(plot_path):
                file_size = os.path.getsize(plot_path)
                ProjectFiles.objects.create(
                    project=project,
                    type='deseq2_visualization',
                    path=plot_path,
                    is_directory=False,
                    file_format='png',
                    size=file_size
                )
                logger.info(f"Registered DESeq2 visualization: {plot_path} with size {file_size} bytes")
                output_files.append(plot_path)
        
        return output_files
    
    except Exception as e:
        logger.error(f"DESeq2 failed: {str(e)}")
        raise RuntimeError(f"DESeq2 failed: {str(e)}")