import argparse
import itertools
import sys

from loguru import logger
import pandas as pd
from pysam import FastaFile
from Bio import Align
from Bio.Align import substitution_matrices


def make_match_only_coords(alignment):
    """Extract the alignment to collect the matched coordinate segments."""
    # This is done by tracing the alignment.path
    # using the code from Bio.Align.PairwiseAlignment.__format__
    # See https://github.com/biopython/biopython/blob/14bf4f2aef448aff47b7b27cf1387eb4272f61db/Bio/Align/__init__.py#L1051  # noqa
    m1_starts = []
    m1_ends = []
    m2_starts = []
    m2_ends = []

    seq1 = alignment.target
    seq2 = alignment.query
    path = alignment.path
    end1, end2 = path[0]

    start1 = end1
    start2 = end2
    for end1, end2 in path[1:]:
        if end1 != start1 and end2 != start2:
            # Match or mismatch
            match_grps = itertools.groupby(
                c1 == c2
                for c1, c2 in zip(seq1[start1:end1], seq2[start2:end2])
            )
            c1_start = start1
            c2_start = start2
            for m, cs in match_grps:
                c_len = len(list(cs))
                c1_end = c1_start + c_len
                c2_end = c2_start + c_len
                # Only take in the matching segments
                if m:
                    m1_starts.append(c1_start)
                    m1_ends.append(c1_end)
                    m2_starts.append(c2_start)
                    m2_ends.append(c2_end)
                c1_start = c1_end
                c2_start = c2_end
        start1 = end1
        start2 = end2
    return m1_starts, m1_ends, m2_starts, m2_ends


def gen_align_coord_str(aligner, source_seq: str, target_seq: str):
    """Align two protein sequences and return the coordinates of the matching segments."""
    alignment = aligner.align(source_seq, target_seq)[0]
    source_starts, source_ends, target_starts, target_ends = make_match_only_coords(alignment)
    return {
        'source_segment_start': source_starts,
        'source_segment_end': source_ends,
        'target_segment_start': target_starts,
        'target_segment_end': target_ends,
    }


def gen_aligner():
    """Create the global sequence aligner

    the parameters are the same as EMBOSS Needle.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Tweak the matrix in BioPython 0.78+ to accomodate amino acid U
    # See https://github.com/biopython/biopython/issues/3205
    # Otherwise use the default matrix and replace U's with X's
    # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    sub_mat = substitution_matrices.load("BLOSUM62")
    sub_mat = sub_mat.select(sub_mat.alphabet + "U")
    # Make U score like X
    sub_mat[:, -1] = -4
    sub_mat[-1, :] = -4
    sub_mat[-1, -1] = 1
    aligner.substitution_matrix = sub_mat
    return aligner


def main(mapping_pth, source_fa_pth, target_fa_pth, out_pth,
         remove_source_id_version=False):
    # Read the mapping table
    logger.info(f"Read mapping table from {mapping_pth}")
    mapping_df = pd.read_table(
        mapping_pth,
        dtype={
            col: 'string'
            for col in [
                'refseq_prot_id', 'hgnc_id', 'source_uniparc_id',
                'uniprot_acc', 'target_uniparc_id', 'mapping_approach'
            ]
        }
    )
    # Select the mappings required sequence alignment
    to_align_df = mapping_df[mapping_df['mapping_approach'] == 'global_seq_align']
    logger.info(f"... total {len(to_align_df):,d} / {len(mapping_df):,d} protein mappings require sequence alignment")

    logger.info(f"Map the sequences from {source_fa_pth} to {target_fa_pth}")
    source_fa = FastaFile(source_fa_pth)
    target_fa = FastaFile(target_fa_pth)

    # Create the global sequence aligner
    aligner = gen_aligner()
    logger.info(f"Global alignment parameters:\n{str(aligner)}")

    # Start alignment
    logger.info("Start global sequence alignment ...")
    aln_results = []
    for i, (source_id, target_id) in enumerate(zip(to_align_df['refseq_prot_id'], to_align_df['uniprot_acc']), 1):
        if remove_source_id_version:
            source_id = source_id[:source_id.rfind('.')]
        # Prior to BioPython 0.78, replace U's since it's an unknown letter to the scoring matrix.
        # source_seq = source_fa.fetch(source_id).replace('U', 'X')
        # target_seq = target_fa.fetch(target_id).replace('U', 'X')
        source_seq = source_fa.fetch(source_id)
        target_seq = target_fa.fetch(target_id)
        coord_str = gen_align_coord_str(aligner, source_seq, target_seq)
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
            ['refseq_prot_id', 'hgnc_id', 'source_uniparc_id',
             'uniprot_acc', 'target_uniparc_id',
             'mapping_approach'],
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
    parser.add_argument("--remove-source-id-version", action='store_true', help="Remove the source entry ID version")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    main(
        args.mapping_tsv,
        args.source_fa,
        args.target_fa,
        args.out_tsv,
        remove_source_id_version=args.remove_source_id_version
    )
