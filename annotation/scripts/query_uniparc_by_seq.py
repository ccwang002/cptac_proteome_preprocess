"""
Get UniParc IDs by the given protein sequence.
"""
import argparse
import asyncio
from itertools import zip_longest
import json
import logging
from pathlib import Path
from typing import List

import aiohttp
import pandas as pd
from pysam import FastaFile
import zstandard as zstd

logger = logging.getLogger(__name__)

_zstd_dctx = zstd.ZstdDecompressor()
_zstd_cctx = zstd.ZstdCompressor(level=9)


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


async def retry_post(session, max_retry=3, **kwargs):
    """Retry POST."""
    for retry in range(max_retry):
        try:
            resp = await session.post(**kwargs)
            resp.raise_for_status()
            return resp
        except aiohttp.ClientResponseError as e:
            logger.error(
                f"Request failed with response {e.status} {e.message} [retry={retry}]"
            )
        except asyncio.TimeoutError:
            logger.error(
                f"Request timed out [retry={retry}]"
            )
        finally:
            if retry == max_retry - 1:
                raise ValueError("Request POST failed")


def parse_uniparc_json(prot_id, j):
    uniparc_id = j["accession"]
    uniparc_checksum = j["sequence"]["checksum"]
    dbrefs = j["dbReference"]
    ref_ids = {}
    for db_type in ["Ensembl", "UniProt", "RefSeq"]:
        ids = set(
            f"{d['id']}.{d['versionI']}"
            for d in dbrefs
            if d["type"].startswith(db_type)
        )
        ref_ids[db_type] = ";".join(ids) or None

    return {
        "original_prot_id": prot_id,
        "uniparc_id": uniparc_id,
        "uniparc_checksum": uniparc_checksum,
        "ensembl_prot_ids": ref_ids["Ensembl"],
        "uniprot_ids": ref_ids["UniProt"],
        "refseq_prot_ids": ref_ids["RefSeq"],
    }


async def query_uniparc(
    prot_id: str, prot_seq: str, session: aiohttp.ClientSession, json_out_pth: Path,
):
    if json_out_pth.exists():
        logger.info(f"... skip {prot_id} query because its JSON exists")
        decompressed = _zstd_dctx.decompress(json_out_pth.read_bytes())
        j = json.loads(decompressed)
        return parse_uniparc_json(prot_id, j)
    try:
        resp = await retry_post(
            session,
            url="https://www.ebi.ac.uk/proteins/api/uniparc/sequence",
            params={
                # 'rfActive': 'true',
                "rfTaxId": "9606",
                "rfDdtype": (
                    "Ensembl,RefSeq,PDB,"
                    "UniProtKB/Swiss-Prot,UniProtKB/Swiss-Prot protein isoforms,"
                    "UniProtKB/TrEMBL"
                ),
            },
            json={"sequence": prot_seq},
        )
    except ValueError:
        logger.error(f"UniParc query of {prot_id} failed")
        return {
            "original_prot_id": prot_id,
            "uniparc_id": None,
            "uniparc_checksum": None,
        }

    j = await resp.json()
    j["sequence"].pop("content")
    j.pop("signatureSequenceMatch", None)
    # Trim the JSON size
    # Write the raw JSON response to external file
    json_out_pth.write_bytes(_zstd_cctx.compress(json.dumps(j).encode()))
    return parse_uniparc_json(prot_id, j)


async def main(
    prot_ids: List[str], protein_fa: FastaFile, out_tsv_pth: Path, json_root: Path
):
    # Make sure all the IDs have their protein sequence in the FASTA
    ids_without_seq = set(prot_ids) - set(protein_fa.references)
    if ids_without_seq:
        logger.critical(
            f'Some IDs are not in the given FASTA: {" ".join(ids_without_seq)}'
        )
        raise ValueError(
            f"Cannot find the sequence of {len(ids_without_seq)} IDs in FASTA"
        )

    logger.info(f"Querying {len(prot_ids)} protein IDs")
    records = []
    ids_without_records = []

    # Group the query by batch
    batch_size = 300
    prot_id_batches = enumerate(grouper(prot_ids, batch_size), 1)
    for batch_ix, prot_ids_one_batch in prot_id_batches:
        # Limit concurrent API calls
        conn = aiohttp.TCPConnector(limit_per_host=10)
        session = aiohttp.ClientSession(
            connector=conn,
            headers={"Accept": "application/json", "Content-type": "application/json"},
        )
        prot_ids_one_batch = [pid for pid in prot_ids_one_batch if pid is not None]
        async with session:
            # Query the IDs
            tasks = []
            for prot_id in prot_ids_one_batch:
                prot_seq = protein_fa.fetch(prot_id)
                json_out_pth = Path(json_dir, f"{prot_id}.json.zst")
                tasks.append(
                    asyncio.create_task(
                        query_uniparc(prot_id, prot_seq, session, json_out_pth),
                        name=prot_id,
                    )
                )

            for i, res in enumerate(asyncio.as_completed(tasks), 1):
                record = await res
                if record['uniparc_id'] is None:
                    ids_without_records.append(record['original_prot_id'])
                records.append(record)
        logger.info(f"... processed {batch_ix * batch_size:,d}/{len(prot_ids):,d} IDs")

    if ids_without_records:
        logger.warning(
            f"These {len(ids_without_records)} IDs have missing query result: "
            f'{" ".join(ids_without_records)}'
        )

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
    parser.add_argument("ids", help="Path to list of IDs (match the FASTA header")
    parser.add_argument("fasta", help="Path to the protein sequence FASTA")
    parser.add_argument("out_tsv", help="Path to the output result TSV")
    parser.add_argument(
        "--json-dir", help="Folder to store the raw query JSONs", default="."
    )

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    with open(args.ids) as f:
        prot_ids = [line.strip() for line in f.read().splitlines() if line]
    protein_fa = FastaFile(args.fasta)
    json_dir = Path(args.json_dir)
    if not json_dir.exists():
        raise SystemExit(f"JSON output folder {json_dir} does not exists.")
    elif not json_dir.is_dir():
        raise SystemExit(f"JSON output folder {json_dir} is not a folder.")
    asyncio.run(main(prot_ids, protein_fa, Path(args.out_tsv), json_dir))
