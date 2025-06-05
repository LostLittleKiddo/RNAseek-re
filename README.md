cd rsa/references/hisat2/featurecounts/

wget ftp://ftp.ensembl.org/pub/release-112/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz
mv Saccharomyces_cerevisiae.R64-1-1.112.gtf Saccharomyces_cerevisiae.R64-1-1.gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz
gunzip gencode.vM35.annotation.gtf.gz
mv gencode.vM35.annotation.gtf.gz Mus_musculus.GRCm39.gtf

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/GRCh38.p14/GFF/ref_GRCh38.p14_top_level.gtf.gz
gunzip ref_GRCh38.p14_top_level.gtf.gz
mv ref_GRCh38.p14_top_level.gtf.gz Homo_sapiens.GRCh38.gtf
