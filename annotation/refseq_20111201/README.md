## RefSeq 20111201
This database has relatively incomplete information.

The genomic information and gene level annotation is missing.


### Workflow
Rename the FASTA header
```bash
sed -E 's/^>gi\|[0-9]+\|ref\|([^|]+)\|/>\1/g' \
    20111201_RefSeq_Human_37-Mouse_37_Trypsin.fasta \
    > 20111201_RefSeq_Human_37-Mouse_37_Trypsin.renamed.fasta
```

Extract the list of protein IDs:

```bash
gzcat 20111201_RefSeq_Human_37-Mouse_37_Trypsin.renamed.fasta.gz \
    | rg '^>.* \[Homo sapiens\]$' \
    | cut -f1 -d' ' \
    | sed 's/^>//g' \
    | sort > refseq_protein_ids.list
```

Run `notebooks/01_make_uniparc_xml_query_tsv.ipynb` to generate the UniParc checksum and download links.

Download the XMLs:

```bash
aria2c -c -j10 --download-result=hide -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    intermediates/refseq_20111201/refseq_protein_ids.list \
    intermediates/refseq_20111201/xmls \
    intermediates/refseq_20111201/refseq_uniparc_mapping.tsv.gz \
    2> intermediates/refseq_20111201/refseq_uniparc_mapping.log
```

Zip the XMLs:
```bash
cd xmls
zstd -T8 -9 *.xml
rm -rf *.xml
cd ..
zip -r -0 xmls.zip xmls/
```