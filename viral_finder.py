import argparse
import gffutils

def main(args):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Viral Finder finds gene neighbourhoods which contain MeldVirus genes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-g", "--gff3", dest="gff_filepath", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="threshold", type=int, required=True, help="Maximum number of gap genes allowed in a neighbourhood.")

    main(parser.parse_args())