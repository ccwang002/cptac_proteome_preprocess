import argparse
import logging
from pathlib import Path
import gzip
import orjson
import xmlschema

logger = logging.getLogger(__name__)


def main(xsd_pth, xml_gz_pth, out_folder):
    logger.info(f'Load XML schema from {xsd_pth}')
    xs = xmlschema.XMLSchema(xsd_pth)

    logger.info(f'Read XML file {xml_gz_pth} ...')
    with gzip.open(xml_gz_pth, 'rt') as f:
        r = xmlschema.XMLResource(f, lazy=True)
        xml_iter = xs.iter_decode(r, path='*')
        for i, entry_d in enumerate(xml_iter, start=1):
            try:
                uniprot_acc = entry_d['accession'][0]
            except TypeError:
                logger.warning(f'Skip an element not an UniProt entry: {entry_d}')
                continue

            out_pth = Path(out_folder, f'{uniprot_acc}.json')
            with out_pth.open('wb') as of:
                of.write(orjson.dumps(entry_d))

            if i % 1000 == 0:
                logger.info(f'... processed {i:,d} entries')

    logger.info(f'Total processed {i:,d} entries')


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
    parser.add_argument("xsd", help="Path to XSD")
    parser.add_argument("xml_gz", help="Path to the gzip'd XML of all UniProt entries")
    parser.add_argument("out_folder", help="Path to the output folder")

    return parser


if __name__ == "__main__":
    parser = setup_cli()
    args = parser.parse_args()
    out_folder = Path(args.out_folder)
    if not out_folder.exists():
        out_folder.mkdir()

    main(args.xsd, args.xml_gz, out_folder)
