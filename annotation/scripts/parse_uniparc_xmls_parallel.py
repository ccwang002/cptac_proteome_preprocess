"""
Get UniParc IDs from the UniParc XMLs.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from loguru import logger
from pathlib import Path
import sys

import pandas as pd
import xmlschema
from xmlschema.etree import ParseError

# This will be initialized in worker processes
xs = None


def initiate_worker(xsd_pth):
    global xs
    xs = xmlschema.XMLSchema(xsd_pth)


def parse_dbref(entry):
    ref_ids = {}
    for db_type in ["Ensembl", "UniProt", "RefSeq"]:
        ids = set()
        for d in entry['dbReference']:
            if d["@type"].startswith(db_type):
                # Skip non-human entries
                d_property = d.get('property', [])
                ncbi_taxid = next(
                    (p['@value'] for p in d_property if p['@type'] == 'NCBI_taxonomy_id'),
                    None
                )
                if ncbi_taxid != '9606':
                    continue

                if '@version' not in d:
                    # Use the UniParc internal version (for UniProt)
                    id_str = f"{d['@id']}.{d['@version_i']}"
                else:
                    id_str = f"{d['@id']}.{d['@version']}"
                ids.add(id_str)
        ref_ids[db_type] = ";".join(ids) or None
    return ref_ids


def process_one_entry(prot_id, xml_pth):
    try:
        data, errors = xs.to_dict(str(xml_pth), validation='lax')
    except ParseError as e:
        logger.error(f"Cannot parse the XML of protein {prot_id}. Likely no record on UniParc")
        # Return an empty record
        record = {
            "original_prot_id": prot_id,
            "uniparc_id": None,
            "uniparc_checksum": None,
            "ensembl_prot_ids": None,
            "uniprot_ids": None,
            "refseq_prot_ids": None,
        }
        return record

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
    return record


def main(
    prot_ids,
    xsd_pth: Path,
    xml_root: Path,
    out_tsv_pth: Path,
    num_processes: int
):

    xml_pths = [Path(xml_root, f'{prot_id}.xml') for prot_id in prot_ids]

    logger.info('Parse XMLs ...')
    records = []

    with ProcessPoolExecutor(max_workers=num_processes, initializer=initiate_worker, initargs=(xsd_pth, )) as executor:
        mapped_jobs = executor.map(
            process_one_entry,
            prot_ids,
            xml_pths,
            chunksize=100
        )
        for i, record_d in enumerate(mapped_jobs, 1):
            records.append(record_d)
            # Show progress
            if i % 1000 == 1:
                logger.info(f"... processed {i:,d} entries")

    # Collect the results
    logger.info(f"Write query results to {out_tsv_pth}")
    result_df = pd.DataFrame.from_records(records)
    result_df.to_csv(out_tsv_pth, sep="\t", index=False)
    logger.info("Done")


def setup_cli():
    logger.remove()  # Remove the default setting

    # Set up the preferred logging colors and format unless overridden by its environment variable
    logger.level("INFO", color="<white>")
    logger.level("DEBUG", color="<d><white>")
    log_format = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> "
        "<b><level>{level: <8}</level></b> "
        "| <level>{message}</level>"
    )
    logger.add(sys.stderr, format=log_format, enqueue=True)

    # Setup CLI parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("ids", help="Path list of IDs")
    parser.add_argument("xsd", help="Path to UniParc XSD")
    parser.add_argument("xml_dir", help="Folder contains the XML files")
    parser.add_argument("out_tsv", help="Path to the output result TSV")
    parser.add_argument("-j", dest="num_processes", type=int, default=1,
                        help="Number of parallel workers")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    with open(args.ids) as f:
        prot_ids = [line.strip() for line in f.read().splitlines() if line]

    main(
        prot_ids,
        args.xsd,
        Path(args.xml_dir),
        Path(args.out_tsv),
        args.num_processes
    )
