# Reference Genomes and Annotations

This directory contains reference genome FASTA files and corresponding gene annotations (GTF/GFF3) for multiple species, suitable for RNA-seq and other genomic workflows.

## Ensembl (Release 114)
- **Danio rerio** (Zebrafish)
  - `Danio_rerio.GRCz11.dna.primary_assembly.fa.gz`: Primary assembly genome
  - `Danio_rerio.GRCz11.114.gtf.gz`: Gene annotations

- **Drosophila melanogaster** (Fruit fly)
  - `Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz`
  - `Drosophila_melanogaster.BDGP6.54.114.gtf.gz`

- **Caenorhabditis elegans** (Nematode)
  - `Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz`
  - `Caenorhabditis_elegans.WBcel235.114.gtf.gz`

- **Mus musculus** (Mouse)
  - `Mus_musculus.GRCm39.dna.toplevel.fa.gz`
  - `Mus_musculus.GRCm39.114.gtf.gz`

- **Homo sapiens** (Human)
  - `Homo_sapiens.GRCh38.dna.toplevel.fa.gz`
  - `Homo_sapiens.GRCh38.114.gtf.gz`

## Ensembl Plants (Release 61)
- **Arabidopsis thaliana**
  - `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`
  - `Arabidopsis_thaliana.TAIR10.61.gff3.gz`

- **Oryza indica** (Rice)
  - `Oryza_indica.ASM465v1.dna.toplevel.fa.gz`
  - `Oryza_indica.ASM465v1.61.gff3.gz`

- **Zea mays** (Maize)
  - `Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz`
  - `Zea_mays.Zm-B73-REFERENCE-NAM-5.0.61.gff3.gz`

## Notes
- GTF files are used for gene structure and splicing information (Ensembl metazoa).
- GFF3 is used for plant genomes (Ensembl Plants).
- FASTA files are used to build indices for tools like HISAT2, STAR, or Salmon.
