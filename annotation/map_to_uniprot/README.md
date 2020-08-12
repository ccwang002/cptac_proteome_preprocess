## Installation

```
conda create -n uniprot python=3.8 notebook flake8 xmlschema zstandard orjson pandas
```

```
fd '.json$' | parallel --bar -j6 'zstd -9 --rm -q {}'
```

```
python 1_extract_xmls.py \
    ../intermediates/uniprot_release_2020_03/uniprot.xsd \
    ../intermediates/uniprot_release_2020_03/proteome_UP000005640_reviewed_yes.xml.gz \
    ../intermediates/uniprot_release_2020_03/parsed_jsons/


python 2_parse_entries.py -j6  \
    ../intermediates/uniprot_release_2020_03/parsed_jsons \
    ../intermediates/uniprot_release_2020_03/uniprot_entries_tsv.gz
```

Make UniParc XML download links:
```python
import pandas as pd

df = pd.read_table('../intermediates/uniprot_release_2020_03/uniprot_entries_tsv.gz')

with open('../intermediates/uniprot_release_2020_03/uniparc_aria2_downloads.links', 'w') as f:
    for t in df.itertuples(index=False):
        print(
            f'https://www.uniprot.org/uniparc/?query=checksum%3A{t.sequence_crc64_checksum}&format=xml\n'
            f'  out=uniparc_xmls/{t.uniprot_acc}.xml',
            file=f
        )
```

Download the UniParc XMLs by:

```bash
aria2c -c -j10 --download-result=hide --console-log-level=warn -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls_parallel.py \
    intermediates/uniprot_release_2020_03/uniprot_entries.list \
    intermediates/uniprot_release_2020_03/uniparc.xsd \
    intermediates/uniprot_release_2020_03/uniparc_xmls \
    intermediates/uniprot_release_2020_03/uniprot_uniparc_mapping.tsv.gz \
    -j 6
```

Compress the XMLs:

```bash
cd xmls
zstd -T8 -9 *.xml
rm -rf *.xml
cd ..
zip -r -0 xmls.zip xmls/
```