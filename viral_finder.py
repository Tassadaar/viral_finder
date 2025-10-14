import sys
import argparse
import gffutils
from Bio import SeqIO

def main(args):
    contigs = list(SeqIO.parse(args.genome_path, "fasta"))
    # we'll read the gff3 in
    db = gffutils.create_db(
        data=args.gff_filepath,
        dbfn=":memory:",
        merge_strategy="create_unique",
        sort_attribute_values=True
    )

    out = []
    neighbourhood = []
    viral_gene_count = 0
    gap_count = 0

    for contig in contigs:

        # THOUGHT: is ordering by start a reasonable assumption?
        for gene in db.region(seqid=contig.id, featuretype="gene"):
            has_name = "Name" in gene.attributes.keys()
            is_viral = "Virus" in str(gene.attributes["Name"]) if has_name else False

            # finding viral genes
            if has_name and is_viral:
                neighbourhood.append(gene)
                viral_gene_count += 1
            elif (has_name and not is_viral) or not has_name:
                gap_count += 1
                if gap_count <= 3:
                    neighbourhood.append(gene)

            # sealing off a neighbourhood
            if gap_count > args.threshold:
                contig = neighbourhood[0].seqid
                start = neighbourhood[0].start
                end = neighbourhood[-1].end
                leng = end - start + 1
                gene_count = len(neighbourhood)

                if viral_gene_count != 0:
                    out.append({
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "length": leng,
                        "gene count": gene_count,
                        "viral gene count": viral_gene_count
                    })

                neighbourhood = []
                viral_gene_count = 0
                gap_count = 0

    for dicty in out:
        print(dicty)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Viral Finder finds gene neighbourhoods which contain MeldVirus genes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-f", "--fasta", dest="genome_path", type=str, required=True)
    parser.add_argument("-g", "--gff3", dest="gff_filepath", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="threshold", type=int, required=True, help="Maximum number of gap genes allowed in a neighbourhood.")

    sys.argv = [
        "viral_finder.py",
        "-f", "/Users/jason/Documents/Roger_Lab/Projects/viral_finder/data/Blastocystis_ST7C_Complete_Assembly.fasta",
        "-g", "/Users/jason/Documents/Roger_Lab/Projects/viral_finder/data/ST7C_gene_predictions_viraltag.gff3",
        "-t", "3"
    ]

    main(parser.parse_args())