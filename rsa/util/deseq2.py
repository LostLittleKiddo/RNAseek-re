# rsa/util/deseq2.py (updated with GSEA, filtering, PDF merging, GFF3 parsing, and 4 significant digits for DESeq2 results)
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
import gseapy as gp
from sklearn.decomposition import PCA
from PyPDF2 import PdfMerger

logger = logging.getLogger(__name__)

def parse_gff3_for_symbols(gff3_path):
    """Parse GFF3 file to extract gene_id to gene_name mapping."""
    try:
        gff3_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        gff3_df = pd.read_csv(gff3_path, sep='\t', comment='#', names=gff3_cols, low_memory=False)
        gff3_df = gff3_df[gff3_df['feature'] == 'gene']
        gene_mapping = {}
        for attr in gff3_df['attributes']:
            gene_id_match = re.search(r'gene_id=([^;]+)', attr)
            name_match = re.search(r'Name=([^;]+)', attr)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                gene_name = name_match.group(1) if name_match else gene_id
                gene_mapping[gene_id] = gene_name
        logger.info(f"Parsed {len(gene_mapping)} gene symbols from GFF3")
        return gene_mapping
    except Exception as e:
        logger.error(f"Failed to parse GFF3 file {gff3_path}: {str(e)}")
        raise RuntimeError(f"Failed to parse GFF3 file: {str(e)}")

def prepare_metadata(meta_data, counts_data):
    """Prepare metadata and align it with counts data."""
    meta_index = meta_data.columns[0]
    meta_data = meta_data.set_index(meta_index)
    shared_samples = counts_data.index.intersection(meta_data.index)
    if len(shared_samples) == 0:
        raise ValueError("No common samples between counts and metadata. Check sample names.")
    logger.info(f"Shared samples: {len(shared_samples)} out of {len(counts_data.index)}")
    counts_data = counts_data.loc[shared_samples]
    meta_data = meta_data.loc[shared_samples]
    return meta_data, counts_data

def create_cluster_heatmap(dds, results_df, output_path, project):
    """Generate a clustered heatmap with dendrograms for significant genes."""
    try:
        sigs = results_df[
            (results_df['padj'] < project.pvalue_cutoff) &
            (results_df['log2FoldChange'].abs() > 1.0) &
            (results_df['baseMean'] > 10.0)
        ]
        if sigs.empty:
            logger.warning("No significant genes found for heatmap")
            raise ValueError("No significant genes to plot in heatmap")
        dds_sigs = dds[:, sigs.index]
        logger.info(f"Subset dds to {len(sigs.index)} significant genes")
        norm_counts = dds_sigs.layers['normed_counts']
        log1p_counts = np.log1p(norm_counts)
        grapher = pd.DataFrame(
            log1p_counts.T,
            index=dds_sigs.var_names,
            columns=dds_sigs.obs_names
        )
        sns.clustermap(
            grapher,
            z_score=0,
            cmap='RdYlBu_r',
            figsize=(10, 8),
            xticklabels=True,
            yticklabels=True,
            cbar_kws={'label': 'Z-score'}
        )
        plt.title(f"Clustered Heatmap - {project.name}")
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

def inspect_deseq2_output(output_file):
    """Inspect DESeq2 results CSV to verify format for GSEA."""
    try:
        df = pd.read_csv(output_file)
        logger.info(f"DESeq2 results shape: {df.shape}")
        logger.info(f"Columns: {df.columns.tolist()}")
        logger.info(f"First few gene IDs: {df.index[:5].tolist()}")
        logger.info(f"Gene symbol column exists: {'gene_symbol' in df.columns}")
        logger.info(f"Number of non-null log2FoldChange: {df['log2FoldChange'].notnull().sum()}")
        logger.info(f"Number of non-null padj: {df['padj'].notnull().sum()}")
    except Exception as e:
        logger.error(f"Failed to inspect DESeq2 output: {str(e)}")
        raise RuntimeError(f"Inspection failed: {str(e)}")

def check_gmt_file(gmt_path):
    """Verify the GMT file format."""
    try:
        with open(gmt_path, 'r') as f:
            lines = f.readlines()
            for line in lines[:5]:
                fields = line.strip().split('\t')
                logger.info(f"GMT line sample: {fields[:3]}")
                if len(fields) < 2:
                    logger.error(f"Invalid GMT format in {gmt_path}")
                    raise RuntimeError(f"Invalid GMT file: {gmt_path}")
        logger.info(f"GMT file appears valid: {gmt_path}")
    except Exception as e:
        logger.error(f"Failed to check GMT file {gmt_path}: {str(e)}")
        raise

def run_gsea(project, deseq2_output_file, gmt_path, output_dir):
    """Run GSEA prerank analysis on DESeq2 results, save filtered outputs, and combine PDFs."""
    try:
        # Read DESeq2 full results
        deseq2_df = pd.read_csv(deseq2_output_file)
        if 'gene_symbol' not in deseq2_df.columns or 'log2FoldChange' not in deseq2_df.columns:
            logger.error("DESeq2 results missing required columns: gene_symbol or log2FoldChange")
            raise RuntimeError("Invalid DESeq2 results format for GSEA")
        
        # Prepare ranked list: gene_symbol and log2FoldChange, convert to uppercase
        ranked_list = deseq2_df[['gene_symbol', 'log2FoldChange']].dropna()
        ranked_list['gene_symbol'] = ranked_list['gene_symbol'].str.upper()  # Convert to uppercase
        ranked_list = ranked_list.set_index('gene_symbol')['log2FoldChange']
        logger.info(f"Prepared ranked list for GSEA with {len(ranked_list)} genes")
        logger.info(f"Sample gene symbols: {ranked_list.index[:5].tolist()}")

        # Run GSEA prerank with adjusted min_size and max_size
        gsea_results = gp.prerank(
            rnk=ranked_list,
            gene_sets=gmt_path,
            outdir=output_dir,
            permutation_num=100,  
            min_size=15,  
            max_size=500, 
            seed=42,
            threads=1
        )
        logger.info(f"GSEA prerank completed. Results saved in: {output_dir}")

        # Filter GSEA results based on specified criteria
        filtered_results = gsea_results.res2d[
            (gsea_results.res2d['NES'].abs() > 1.5) &
            (gsea_results.res2d['FDR q-val'] < 0.25) &
            (gsea_results.res2d['NOM p-val'] < 0.05)
        ]
        logger.info(f"Filtered GSEA results to {len(filtered_results)} gene sets with "
                    f"|NES| > 1.5, FDR < 0.25, NOM p-val < 0.05")

        # Drop 'Name' column and save filtered GSEA results to CSV without index
        filtered_results = filtered_results.drop(columns=['Name'], errors='ignore')
        gsea_output_file = os.path.join(output_dir, "gsea_results.csv")
        filtered_results.to_csv(gsea_output_file, index=False)
        logger.info(f"Filtered GSEA results CSV saved to: {gsea_output_file}")

        # Combine PDFs in output_dir/prerank into gsea_plot.pdf
        prerank_dir = os.path.join(output_dir, "prerank")
        combined_pdf_path = os.path.join(output_dir, "gsea_plot.pdf")
        merger = PdfMerger()
        pdf_files = [f for f in os.listdir(prerank_dir) if f.endswith('.pdf')] if os.path.exists(prerank_dir) else []
        
        if pdf_files:
            for pdf_file in sorted(pdf_files):  # Sort for consistent order
                pdf_path = os.path.join(prerank_dir, pdf_file)
                try:
                    merger.append(pdf_path)
                    logger.info(f"Added {pdf_path} to combined PDF")
                except Exception as e:
                    logger.warning(f"Failed to append {pdf_path}: {str(e)}")
            try:
                merger.write(combined_pdf_path)
                logger.info(f"Combined GSEA PDFs saved to: {combined_pdf_path}")
            except Exception as e:
                logger.error(f"Failed to save combined PDF: {str(e)}")
                combined_pdf_path = None
            finally:
                merger.close()
        else:
            logger.warning(f"No PDF files found in {prerank_dir}")
            combined_pdf_path = None

        # Register GSEA outputs in ProjectFiles
        output_files = [gsea_output_file]
        
        if combined_pdf_path and os.path.exists(combined_pdf_path):
            file_size = os.path.getsize(combined_pdf_path)
            ProjectFiles.objects.create(
                project=project,
                type='gsea_visualization',
                path=combined_pdf_path,
                is_directory=False,
                file_format='pdf',
                size=file_size
            )
            logger.info(f"Registered combined GSEA PDF: {combined_pdf_path} with size {file_size} bytes")
            output_files.append(combined_pdf_path)

        for file_path in [gsea_output_file]:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
                file_type = 'gsea_output'
                file_format = 'csv'
                ProjectFiles.objects.create(
                    project=project,
                    type=file_type,
                    path=file_path,
                    is_directory=False,
                    file_format=file_format,
                    size=file_size
                )
                logger.info(f"Registered GSEA file: {file_path} with size {file_size} bytes")

        return output_files

    except Exception as e:
        logger.error(f"GSEA failed: {str(e)}")
        # Return empty list to allow DESeq2 to complete
        return []

def run_deseq2(project, counts_file, metadata_file, output_dir):
    """
    Run DESeq2 on counts.csv and metadata.csv, adding gene symbols from GFF3.
    Filter results by project.pvalue_cutoff, log2FoldChange > 1, and baseMean > 10.
    Generate cluster heatmap, PCA, and GSEA analysis.

    Args:
        project: Project instance (contains species and pvalue_cutoff).
        counts_file: Path to counts.csv from FeatureCounts.
        metadata_file: Path to metadata.csv.
        output_dir: Directory for DESeq2 and GSEA output (deseq2_results.csv, gsea_results.csv, and plots).

    Returns:
        list: Paths to deseq2_results.csv, gsea_results.csv, and visualization PNGs/PDFs.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "deseq2_results.csv")
    full_output_file = os.path.join(output_dir, "deseq2_full_results.csv")
    heatmap_output = os.path.join(output_dir, "heatmap.png")
    pca_output = os.path.join(output_dir, "pca_plot.png")

    species_to_gff3 = {
        'human': 'Homo_sapiens.GRCh38.gff3',
        'mouse': 'Mus_musculus.GRCm39.gff3',
        'yeast': 'Saccharomyces_cerevisiae.R64-1-1.gff3'
    }
    gff3_file = species_to_gff3.get(project.species.lower(), None)
    if not gff3_file:
        raise RuntimeError(f"No GFF3 file defined for species: {project.species}")
    gff3_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'gff3', gff3_file)
    if not os.path.exists(gff3_path):
        raise RuntimeError(f"GFF3 file not found: {gff3_path}")

    try:
        counts = pd.read_csv(counts_file, sep='\t')
        counts = counts.set_index('Geneid')
        columns_to_remove = ['Chr', 'Start', 'End', 'Strand', 'Length']
        counts = counts.drop(columns=columns_to_remove, errors='ignore')
        numeric_counts = counts.select_dtypes(include=['int64', 'float64'])
        counts = counts[numeric_counts.sum(axis=1) > 0]
        counts = counts.T

        metadata = pd.read_csv(metadata_file)
        metadata, counts = prepare_metadata(metadata, counts)
        
        dds = DeseqDataSet(counts=counts, metadata=metadata, design='condition')
        dds.deseq2()
        stat_res = DeseqStats(dds, n_cpus=1, contrast=['condition', metadata['condition'].unique()[0], metadata['condition'].unique()[1]])
        stat_res.summary()
        results_df = stat_res.results_df

        gene_mapping = parse_gff3_for_symbols(gff3_path)
        results_df['gene_symbol'] = results_df.index.map(gene_mapping)
        
        # Round numerical columns to 4 significant digits for full results
        numeric_cols = results_df.select_dtypes(include=['float64', 'int64']).columns
        results_df[numeric_cols] = results_df[numeric_cols].round(4)
        results_df.to_csv(full_output_file)
        logger.info(f"Full DESeq2 results saved to: {full_output_file}")
        if os.path.exists(full_output_file):
            file_size = os.path.getsize(full_output_file)
            ProjectFiles.objects.create(
                project=project,
                type='deseq2_full_output',
                path=full_output_file,
                is_directory=False,
                file_format='csv',
                size=file_size
            )
            logger.info(f"Registered full DESeq2 output CSV: {full_output_file} with size {file_size} bytes")

        pvalue_cutoff = project.pvalue_cutoff
        results_df = results_df[
            (results_df['padj'] < pvalue_cutoff) &
            (results_df['log2FoldChange'].abs() > 1.0) &
            (results_df['baseMean'] > 10.0)
        ]
        logger.info(f"Filtered DESeq2 results to {len(results_df)} genes with padj < {pvalue_cutoff}, "
                    f"|log2FoldChange| > 1.0, and baseMean > 10.0")
        
        # Round numerical columns to 4 significant digits for filtered results
        numeric_cols = results_df.select_dtypes(include=['float64', 'int64']).columns
        results_df[numeric_cols] = results_df[numeric_cols].round(4)
        results_df.to_csv(output_file)
        logger.info(f"DESeq2 results saved to: {output_file}")

        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            ProjectFiles.objects.create(
                project=project,
                type='deseq_output',
                path=output_file,
                is_directory=False,
                file_format='csv',
                size=file_size
            )
            logger.info(f"Registered DESeq2 output CSV: {output_file} with size {file_size} bytes")

        output_files = [output_file]
        create_pca_plot(counts, metadata, pca_output, project)
        create_cluster_heatmap(dds, results_df, heatmap_output, project)

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
        
        inspect_deseq2_output(output_file)

        # species_to_gmt = {
        #     'human': 'c5.go.v2023.1.Hs.symbols.gmt',
        #     'mouse': 'm2.all.v2023.1.Mm.symbols.gmt',
        #     'yeast': 'c5.go.bp.v2023.1.Sc.symbols.gmt'
        # }
        # gmt_file = species_to_gmt.get(project.species.lower())
        # if gmt_file:
        #     gmt_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'gsea', gmt_file)
        #     if os.path.exists(gmt_path):
        #         check_gmt_file(gmt_path)
        #         # Use full DESeq2 results for GSEA
        #         gsea_output_files = run_gsea(project, full_output_file, gmt_path, output_dir)
        #         output_files.extend(gsea_output_files)
        #     else:
        #         logger.warning(f"GMT file not found: {gmt_path}")
        # else:
        #     logger.warning(f"No GMT file defined for species: {project.species}")

        return output_files
    
    except Exception as e:
        logger.error(f"DESeq2 failed: {str(e)}")
        raise RuntimeError(f"DESeq2 failed: {str(e)}")