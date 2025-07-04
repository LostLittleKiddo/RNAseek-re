#!/bin/bash
set -e

mkdir -p rsa/references/{gtf,index}

# ---------- GTF/GFF3 FILES ----------
cd rsa/references/gtf

# Zebrafish (Danio rerio)
wget -c https://ftp.ensembl.org/pub/release-114/gtf/danio_rerio/Danio_rerio.GRCz11.114.gtf.gz

# Drosophila melanogaster
wget -c https://ftp.ensembl.org/pub/release-114/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.114.gtf.gz

# Caenorhabditis elegans
wget -c https://ftp.ensembl.org/pub/release-114/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.114.gtf.gz

# Mus musculus (Mouse)
wget -c https://ftp.ensembl.org/pub/release-114/gtf/mus_musculus/Mus_musculus.GRCm39.114.gtf.gz

# Homo sapiens (Human)
wget -c https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

# Arabidopsis thaliana
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gff3.gz

# Oryza sativa
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.61.gff3.gz

# Zea mays
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.61.gff3.gz

# Saccharomyces cerevisiae (Yeast)
wget -c https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz

cd ../index

# ---------- FASTA & HISAT2 INDEX ----------
declare -A species=(
  [zebrafish]="https://ftp.ensembl.org/pub/release-114/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
  [fly]="https://ftp.ensembl.org/pub/release-114/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz"
  [worm]="https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
  [mouse]="https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/dna_index/Mus_musculus.GRCm39.dna.toplevel.fa.gz"
  [human]="https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
  [arabidopsis]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/arabidopsis_thaliana/dna_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
  [oryza]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/oryza_sativa/dna_index/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
  [maize]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/zea_mays/dna_index/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
  [yeast]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-61/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
)

for name in "${!species[@]}"; do
  mkdir -p "$name"
  cd "$name"

  fasta_url="${species[$name]}"
  fasta_file=$(basename "$fasta_url")
  base_name="${fasta_file%.gz}"

  echo "📥 Downloading $name genome..."
  wget -c "$fasta_url"

  echo "📦 Unzipping $fasta_file..."
  gunzip -f "$fasta_file"

  echo "🔧 Building HISAT2 index for $base_name..."
  hisat2-build "$base_name" "${name}_hisat2_index"

  cd ..
done

echo "✅ All reference genomes downloaded and indexed."
