## RefSeq 20160914

Note that mitochondrial genes (`YP_*`) are reported on chrMT ([rCRS] sequence; 16,569bp) instead of chrM (official GRCh37; 16,571bp). chrMT is identical to hg38's chrM.

[rCRS]: https://en.wikipedia.org/wiki/Cambridge_Reference_Sequence


### Workflow
Use `notebooks/01_cptac_refseq_granges_annotation.Rmd` to obtain the full RefSeq protein ID list.

Run `notebooks/02_make_uniparc_xml_query_tsv.ipynb` to generate the CRC64 checksum for the missing RefSeq IDs.

Use the checksum to query UniProt to get the UniParc XMLs:

```bash
aria2c -c -j10 --download-result=hide --console-log-level=warn -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    intermediates/refseq_20160914/refseq_protein_ids.list \
    intermediates/refseq_20160914/xmls \
    intermediates/refseq_20160914/refseq_uniparc_mapping.tsv.gz \
    2> intermediates/refseq_20160914/refseq_uniparc_mapping.log
```

Compress the XMLs:

```bash
cd xmls
zstd -T8 -9 *.xml
rm -rf *.xml
cd ..
zip -r -0 xmls.zip xmls/
```

Run `notebooks/01_cptac3_refseq_granges_annotation.Rmd` to merge the RefSeq annotations with UniParc.