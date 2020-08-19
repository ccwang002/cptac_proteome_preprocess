import argparse
from concurrent.futures import ProcessPoolExecutor
import sys
from loguru import logger
from pathlib import Path
import pandas as pd
import orjson
import zstandard as zstd


# This will be initialized in worker processes
_zstd_dctx = None


def initiate_worker():
    global _zstd_dctx
    _zstd_dctx = zstd.ZstdDecompressor()


def read_and_parse_one_entry(json_pth):
    entry_d = orjson.loads(_zstd_dctx.decompress(
        open(json_pth, 'rb').read()
    ))
    try:
        out_d = parse_uniprot_entry(entry_d)
        return out_d
    except Exception:
        logger.opt(exception=True).error(f'Failed to parse {json_pth.name}')


def parse_uniprot_entry(entry_d):
    uniprot_acc = entry_d['accession'][0]
    uniprot_name = entry_d['name'][0]

    # parse protein name
    protein_name = entry_d['protein']['recommendedName']['fullName']
    if not isinstance(protein_name, str):
        protein_name = protein_name['$']

    # parse gene name
    try:
        gene_name = entry_d['gene'][0]['name'][0]['$']
    except KeyError:
        gene_name = None

    # parse HGNC
    # note that one UniProt entry can have multiple HGNC mappings
    # e.g. Q9BQY6
    hgnc_ds = [d for d in entry_d['dbReference'] if d['@type'] == 'HGNC']
    if hgnc_ds:
        if len(hgnc_ds) == 1:
            hgnc_d = hgnc_ds[0]
            hgnc_gene_id = hgnc_d['@id']
            hgnc_gene_name = next(
                (d['@value'] for d in hgnc_d['property'] if d['@type'] == 'gene designation'),
                None
            )
        else:
            hgnc_gene_id = ";".join([d['@id'] for d in hgnc_ds])
            hgnc_gene_name = ";".join([
                next(
                    (d['@value'] for d in hgnc_d['property'] if d['@type'] == 'gene designation'),
                    ""  # keep the empty value
                )
                for hgnc_d in hgnc_ds
            ])
            logger.warning(f"Entry {uniprot_acc} has multiple HGNC mappings: {hgnc_gene_id}")
    else:
        hgnc_gene_id = None
        hgnc_gene_name = None

    seq_d = entry_d['sequence']
    return {
        'uniprot_acc': uniprot_acc,
        'uniprot_name': uniprot_name,
        'uniprot_protein_name': protein_name,
        'uniprot_gene_name': gene_name,
        'hgnc_gene_ids': hgnc_gene_id,
        'hgnc_gene_names': hgnc_gene_name,
        'sequence_modified_date': seq_d['@modified'],
        'sequence_crc64_checksum': seq_d['@checksum'],
        'sequence': seq_d['$']
    }


def main(json_folder: Path, out_tsv_pth: Path, num_processes: int):
    logger.info(f'Load JSONs from {json_folder}')
    json_pths = sorted(json_folder.glob('*.json.zst'))
    logger.info(f'... found {len(json_pths):,d} entries to process')

    logger.info('Parse JSONs ...')
    results = []

    with ProcessPoolExecutor(max_workers=num_processes, initializer=initiate_worker) as executor:
        mapped_jobs = executor.map(
            read_and_parse_one_entry,
            json_pths,
            chunksize=100
        )
        for i, result_d in enumerate(mapped_jobs, 1):
            results.append(result_d)
            if i % 1000 == 0:
                logger.info(f'... processed {i:6,d}/{len(json_pths):,d} entries')
    logger.info(f'Parsed all JSONs. Write result to {out_tsv_pth}')

    df = pd.DataFrame.from_records(
        results,
        columns=[
            'uniprot_acc', 'uniprot_name',
            'uniprot_protein_name', 'uniprot_gene_name',
            'hgnc_gene_ids', 'hgnc_gene_names',
            'sequence_modified_date', 'sequence_crc64_checksum',
            'sequence'
        ]
    )
    df.to_csv(out_tsv_pth, index=False, sep='\t')
    logger.info('Done')


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
    parser.add_argument("json_folder", help="Path to XSD")
    parser.add_argument("out_tsv", help="Path to the output TSV")
    parser.add_argument("-j", dest="num_processes", type=int, default=1,
                        help="Number of parallel workers")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    main(
        Path(args.json_folder),
        Path(args.out_tsv),
        args.num_processes
    )
