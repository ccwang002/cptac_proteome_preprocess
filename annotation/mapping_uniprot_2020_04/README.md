## UniProt release 2020_04 reviewed human proteome
The entries in use are defined by this query:

    reviewed:yes AND organism:"Homo sapiens (Human) [9606]" AND proteome:up000005640


### Download UniProtKB XMLs
Download `proteome_UP000005640_reviewed_yes.xml.gz` and `proteome_UP000005640_reviewed_yes.fasta.gz`
using the UniProtKB web interface:

```bash
curl -Lo proteome_UP000005640_reviewed_yes.fasta.gz \
    'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+proteome%3Aup000005640&sort=score&format=fasta&compress=yes'

curl -Lo proteome_UP000005640_reviewed_yes.xml.gz \
    'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+proteome%3Aup000005640&sort=score&format=xml&compress=yes'

curl -LO https://www.uniprot.org/docs/uniparc.xsd
curl -LO https://www.uniprot.org/docs/uniprot.xsd
```


Parse the giant UniProtKB XML:

```bash
# Split the giant XML into one JSON per entry
python 1_extract_xmls.py \
    ../intermediates/uniprot_release_2020_04/uniprot.xsd \
    ../intermediates/uniprot_release_2020_04/proteome_UP000005640_reviewed_yes.xml.gz \
    ../intermediates/uniprot_release_2020_04/parsed_jsons/

# Compress the JSONs
cd ../intermediates/uniprot_release_2020_04/parsed_jsons/
fd '.json$' . | parallel --bar -j8 'zstd -9 --rm -q {}'

# Parse the JSONs
python 2_parse_entries.py -j8  \
    ../intermediates/uniprot_release_2020_04/parsed_jsons \
    ../intermediates/uniprot_release_2020_04/uniprot_entries_tsv.gz

# Package all JSONs
cd ../intermediates/uniprot_release_2020_04/
zip -r -0 parsed_jsons.zip parsed_jsons
rm -rf parsed_jsons
```


### Add UniParc annotation
Generate UniParc XML download links:

```python
import pandas as pd

df = pd.read_table('../intermediates/uniprot_release_2020_04/uniprot_entries_tsv.gz')

with open('../intermediates/uniprot_release_2020_04/uniparc_aria2_downloads.links', 'w') as f:
    for t in df.itertuples(index=False):
        # Alternatively, use this URL:
        # 'https://www.uniprot.org/uniparc/?query=uniprot%3A{t.uniprot_acc}&sort=score&direct=yes&format=xml'
        print(
            f'https://www.uniprot.org/uniparc/?query=checksum%3A{t.sequence_crc64_checksum}&format=xml\n'
            f'  out=uniparc_xmls/{t.uniprot_acc}.xml',
            file=f
        )
```

Download the UniParc XMLs by:

```bash
cd ../intermediates/uniprot_release_2020_04
aria2c -c -j10 --retry-wait=10 --download-result=hide --console-log-level=warn -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
gzip -d -c uniprot_entries_tsv.gz | tail -n+2 | cut -f1 > uniprot_entries.list

python scripts/parse_uniparc_xmls_parallel.py \
    intermediates/uniprot_release_2020_04/uniprot_entries.list \
    intermediates/uniprot_release_2020_04/uniparc.xsd \
    intermediates/uniprot_release_2020_04/uniparc_xmls \
    intermediates/uniprot_release_2020_04/uniprot_uniparc_mapping.tsv.gz \
    -j 8
```

Compress the XMLs:

```bash
cd uniparc_xmls
fd '.xml$' | parallel --bar -j6 'zstd -9 -q  --rm {}'
cd ..
zip -r -0 uniparc_xmls.zip uniparc_xmls/
```

Run `notebooks/01_combine_annotations.Rmd`.



## Create mapping to RefSeq annotations
Run `notebooks/{02-05}_map_refseq*.Rmd` to create the mappings.

Run sequence alignment per RefSeq annotation version:

```bash
# RefSeq 2018
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20180629_to_uniprot_2020_04_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.0/DCC/RefSeq_20180629/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2020_04.fasta.gz  \
    tracked_results/mappings/refseq_20180629_to_uniprot_2020_04_mappable.coord_mapping.tsv.gz


# RefSeq 2016
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20160914_to_uniprot_2020_04_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.0/DCC/RefSeq_20160914/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2020_04.fasta.gz  \
    tracked_results/mappings/refseq_20160914_to_uniprot_2020_04_mappable.coord_mapping.tsv.gz


# RefSeq 2013
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20130727_to_uniprot_2020_04_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.0/DCC/RefSeq_20130727/RefSeq.20130727-Human.contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2020_04.fasta.gz  \
    tracked_results/mappings/refseq_20130727_to_uniprot_2020_04_mappable.coord_mapping.tsv.gz \
    --remove-source-id-version


# RefSeq 2011
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20111201_to_uniprot_2020_04_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.0/DCC/RefSeq_20111201/20111201_RefSeq_Human_37-Mouse_37_Trypsin.renamed.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2020_04.fasta.gz  \
    tracked_results/mappings/refseq_20111201_to_uniprot_2020_04_mappable.coord_mapping.tsv.gz
```
