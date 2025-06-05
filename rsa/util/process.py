# rsa/util/preprocess_gff.py
import re
import logging
import os

logger = logging.getLogger(__name__)

def preprocess_gff(input_gff, output_gff):
    """
    Add gene attribute to exon features based on Parent attribute in GFF3 file.
    Use Parent ID as fallback for non-coding RNAs lacking gene/Name attributes.
    
    Args:
        input_gff: Path to input GFF file.
        output_gff: Path to output preprocessed GFF file.
    """
    # Map parent IDs to gene names or IDs
    parent_id_map = {}
    with open(input_gff, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            if cols[2] in ('gene', 'mRNA', 'tRNA', 'rRNA', 'ncRNA'):
                attrs = dict(re.findall(r'(\w+)=([^;]+)', cols[8]))
                gene_name = attrs.get('gene', attrs.get('Name', attrs.get('ID', '')))
                parent_id_map[attrs.get('ID', '')] = gene_name or attrs.get('ID', '')

    # Rewrite GFF file, adding gene attribute to exon features
    with open(input_gff, 'r') as f, open(output_gff, 'w') as out:
        for line in f:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            cols = line.strip().split('\t')
            if cols[2] == 'exon':
                attrs = dict(re.findall(r'(\w+)=([^;]+)', cols[8]))
                parent_id = attrs.get('Parent', '')
                if parent_id in parent_id_map:
                    attrs['gene'] = parent_id_map[parent_id]
                    cols[8] = ';'.join(f"{k}={v}" for k, v in attrs.items())
                else:
                    logger.warning(f"No parent found for exon with Parent={parent_id}")
            out.write('\t'.join(cols) + '\n')
    logger.info(f"Preprocessed GFF file written to: {output_gff}")

if __name__ == "__main__":
    input_gff = "/home/littlekiddo/Desktop/RNAseek-re/rsa/util/annotation/Saccharomyces_cerevisiae.R64-1-1.gff"
    output_gff = "/home/littlekiddo/Desktop/RNAseek-re/rsa/util/annotation/Saccharomyces_cerevisiae.R64-1-1.processed.gff"
    preprocess_gff(input_gff, output_gff)