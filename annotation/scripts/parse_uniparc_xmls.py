"""
Get UniParc IDs from the UniParc XMLs.
"""
import argparse
import logging
from pathlib import Path

import pandas as pd
import xmlschema

logger = logging.getLogger(__name__)


def parse_dbref(entry):
    ref_ids = {}
    for db_type in ["Ensembl", "UniProt", "RefSeq"]:
        ids = set()
        for d in entry['dbReference']:
            if d["@type"].startswith(db_type):
                if '@version' not in d:
                    # Use the internal version
                    # logger.warning(f'This dbRef has no external version: {d!r}')
                    id_str = f"{d['@id']}.{d['@version_i']}"
                else:
                    id_str = f"{d['@id']}.{d['@version']}"
                ids.add(id_str)
        ref_ids[db_type] = ";".join(ids) or None
    return ref_ids


def main(
    prot_ids,
    xml_root: Path,
    out_tsv_pth: Path
):
    xs = xmlschema.XMLSchema('https://www.uniprot.org/docs/uniparc.xsd')

    records = []
    for prot_id in prot_ids:
        xml_pth = Path(xml_root, f'{prot_id}.xml')
        data, errors = xs.to_dict(str(xml_pth), validation='lax')
        if len(data['entry']) > 1:
            ids = [x['accession'] for x in data['entry']]
            logger.warning(f'{prot_id} got multiple uniparc entries: {" ".join(ids)}')

        entry = data['entry'][0]
        uniparc_id = entry['accession']
        uniparc_cksum = entry['sequence']['@checksum']
        ref_ids = parse_dbref(entry)
        record = {
            "original_prot_id": prot_id,
            "uniparc_id": uniparc_id,
            "uniparc_checksum": uniparc_cksum,
            "ensembl_prot_ids": ref_ids["Ensembl"],
            "uniprot_ids": ref_ids["UniProt"],
            "refseq_prot_ids": ref_ids["RefSeq"],
        }
        records.append(record)

    # Collect the results
    logger.info(f"Write query results to {out_tsv_pth}")
    result_df = pd.DataFrame.from_records(records)
    result_df.to_csv(out_tsv_pth, sep="\t", index=False)


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = "[%(asctime)s][%(levelname)-7s] %(message)s"
    log_formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
    console.setFormatter(log_formatter)

    # Setup CLI parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("ids", help="Path list of IDs")
    parser.add_argument("xml_dir", help="Folder contains the XML files")
    parser.add_argument("out_tsv", help="Path to the output result TSV")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    with open(args.ids) as f:
        prot_ids = [line.strip() for line in f.read().splitlines() if line]

    main(
        prot_ids,
        Path(args.xml_dir),
        Path(args.out_tsv)
    )
