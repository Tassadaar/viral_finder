import sys
import argparse
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def main(args):
    contigs: list[SeqRecord] = list(SeqIO.parse(args.genome_path, "fasta"))
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
    # gene pattern (p=penton, b=polB, h=hypothetical, F=FtsK, a=adenain, e=putative fiber e, i=integrase, r=non-viral labeled genes)
    pattern_lookup = {
        "penton": "p",
        "PolB": "b",
        "hypothetical": "h",
        "FtsK": "f",
        "adenain": "a",
        "putative": "e",
        "integrase": "i",
    }
    pattern = []

    for contig in contigs:
        print(f"Processing {contig.id}...\n")

        # THOUGHT: is ordering by start a reasonable assumption?
        genes = list(db.region(seqid=contig.id, featuretype="gene", completely_within=True))

        for gene in genes:
            has_name = "Name" in gene.attributes.keys()
            is_viral = "Virus" in str(gene.attributes["Name"]) if has_name else False

            # skip non-viral starts
            if (not has_name or (has_name and not is_viral)) and len(neighbourhood) == 0:
                continue

            # finding viral genes
            if has_name and is_viral:
                neighbourhood.append(gene)
                viral_gene_count += 1

                # resetting gap count if viral gene is reach before threshold
                if gap_count == args.threshold: gap_count = 0

                # adding pattern
                for key, value in pattern_lookup.items():
                    if key in str(gene.attributes["Name"]):
                        pattern.append(value)

            # finding non-viral genes in gaps
            elif (has_name and not is_viral) or not has_name:
                gap_count += 1
                if gap_count <= 3:
                    neighbourhood.append(gene)
                    pattern.append("r")

            # sealing off a neighbourhood when threshold is reached or at the end of a contig
            if gap_count > args.threshold or gene == genes[-1]:
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
                        "viral gene count": viral_gene_count,
                        "pattern": "".join(pattern)
                    })

                neighbourhood = []
                viral_gene_count = 0
                gap_count = 0
                pattern = []

    pd.DataFrame(out).to_csv("result.csv", index=False)
    print("Results have been saved to ./result.csv.")
    print("Exiting...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Viral Finder finds gene neighbourhoods which contain MeldVirus genes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-f", "--fasta", dest="genome_path", type=str, required=True)
    parser.add_argument("-g", "--gff3", dest="gff_filepath", type=str, required=True)
    parser.add_argument("-t", "--threshold", dest="threshold", type=int, required=True, help="Maximum number of gap genes allowed in a neighbourhood.")

    main(parser.parse_args())