## GENCODE 34 (corresponding to Ensembl 100)

### Workflow
Download the source files:

```bash
mkdir intermediates/gencode_34

# GENCODE protein sequence
curl -L http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_translations.fa.gz -o intermediates/gencode_34/gencode.v34.pc_translations.fa.gz

curl -L https://www.uniprot.org/docs/uniparc.xsd -o intermediates/gencode_34/uniparc.xsd

# EnsDb (Ensembl 100) from AnnotationHub AH79689
curl -L http://s3.amazonaws.com/annotationhub/AHEnsDbs/v100/EnsDb.Hsapiens.v100.sqlite -o ~/Box/Resources/EnsDb.Hsapiens.v100.sqlite
```

Run `notesbooks/01_merge_same_protein_seq.Rmd`

```bash
# Compress and index the fasta file
bgzip -l9 gencode_34_unique_protein_seq.fasta
samtools faidx gencode_34_unique_protein_seq.fasta.gz
```

Run `notebooks/02_gen_uniparc_checksums.ipynb`

Download the uniparc XMLs:

```bash
cd ../intermediates/gencode_34
aria2c -c -j10 \
    --auto-file-renaming=false \
    --download-result=hide \
    --console-log-level=warn \
    -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls_parallel.py \
    intermediates/gencode_34/gencode_34_unique_protein_entries.list \
    intermediates/gencode_34/uniparc.xsd \
    intermediates/gencode_34/xmls \
    intermediates/gencode_34/gencode_34_uniparc_mapping.tsv.gz \
    -j 8
```

