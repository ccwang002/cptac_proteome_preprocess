import argparse
import importlib
import itertools
import sys

from loguru import logger
import pandas as pd
from pysam import FastaFile
from Bio import Align
from Bio.Align import substitution_matrices

m = importlib.import_module('3_run_alignment')

def main(mapping_pth, source_fa_pth, target_fa_pth, out_pth):
    mapping_cols = [
        'ensembl_prot_id', 'uniq_protein_seq_entry_no',
        'symbol', 'hgnc_id',
        'source_uniparc_id',
        'uniprot_acc', 'target_uniparc_id', 'mapping_approach'
    ]
    # Read the mapping table
    logger.info(f"Read mapping table from {mapping_pth}")
    mapping_df = pd.read_table(
        mapping_pth,
        dtype={
            col: 'string'
            for col in mapping_cols
        }
    )
    # Select the mappings required sequence alignment
    to_align_df = mapping_df[mapping_df['mapping_approach'] == 'global_seq_align']
    logger.info(f"... total {len(to_align_df):,d} / {len(mapping_df):,d} protein mappings require sequence alignment")

    logger.info(f"Map the sequences from {source_fa_pth} to {target_fa_pth}")
    source_fa = FastaFile(source_fa_pth)
    target_fa = FastaFile(target_fa_pth)

    # Create the global sequence aligner
    aligner = m.gen_aligner()
    logger.info(f"Global alignment parameters:\n{str(aligner)}")

    # Start alignment
    logger.info("Start global sequence alignment ...")
    aln_results = []
    for i, (source_id, target_id) in enumerate(zip(to_align_df['uniq_protein_seq_entry_no'], to_align_df['uniprot_acc']), 1):
        # Prior to BioPython 0.78, replace U's since it's an unknown letter to the scoring matrix.
        # source_seq = source_fa.fetch(source_id).replace('U', 'X')
        # target_seq = target_fa.fetch(target_id).replace('U', 'X')
        source_seq = source_fa.fetch(source_id)
        target_seq = target_fa.fetch(target_id)
        coord_str = m.gen_align_coord_str(aligner, source_seq, target_seq)
        aln_results.append(coord_str)
        if i % 1000 == 0:
            logger.info(f"... processed {i:d} / {len(to_align_df)} mappings")
    logger.info("Global sequence alignment complete")

    # Collect all alignments together
    logger.info("Organize the results ...")
    align_coord_series = pd.DataFrame.from_records(aln_results, index=to_align_df.index)
    out_df = (
        mapping_df
        # Not all the mappings require alignment, so join by the index here
        .join(align_coord_series)
        # Flatten the nested coordinate segments
        .set_index(
            mapping_cols,
            verify_integrity=True
        )
        .apply(pd.Series.explode)
        .apply(pd.Series.astype, dtype=pd.Int64Dtype())
        .reset_index()
    )
    # Convert the coordinates to be 1-based close interval
    out_df['source_segment_start'] += 1
    out_df['target_segment_start'] += 1
    # Calculate the coordinate change distance for fast conversion
    out_df['target_to_source_distance'] = out_df['target_segment_start'] - out_df['source_segment_start']

    logger.info(f"Export the results to {out_pth}")
    out_df.to_csv(out_pth, sep='\t', index=False)
    logger.success("Done")


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
    logger.add(sys.stderr, format=log_format)

    # Setup CLI parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("mapping_tsv", help="TSV table of the protein ID mapping")
    parser.add_argument("source_fa", help="Path to the source protein FASTA")
    parser.add_argument("target_fa", help="Path to the target protein FASTA")
    parser.add_argument("out_tsv", help="Output coordinate mapping TSV")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    main(
        args.mapping_tsv,
        args.source_fa,
        args.target_fa,
        args.out_tsv,
    )
