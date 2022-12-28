#!/bin/env python3
import shutil, copy
from Bio import SeqIO
from BCBio import GFF

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--reference', help='Genbank accession of reference sequence (will be fetched from genbank)')
    parser.add_argument('--output-dir', required=True, type=str, help='Output directory')
    return parser.parse_args()

def get_reference_sequence(accession):
    from Bio import Entrez
    Entrez.email = "hello@nextstrain.org"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    print(f"Fetching reference sequence {accession} from genbank")
    return SeqIO.read(handle, "genbank")


if __name__=="__main__":
    args = parse_args()
    reference = get_reference_sequence(args.reference)

    # copy the template to the output directory
    shutil.copytree(f"files", f"{args.output_dir}/files")

    # write the reference sequence to the output directory
    SeqIO.write(reference, f"{args.output_dir}/files/reference.fasta", "fasta")

    new_features = []
    for feature in reference.features:
        if feature.type == "CDS":
            new_feature = copy.deepcopy(feature)
            new_feature.type = "gene"
            new_feature.qualifiers["gene_name"] = new_feature.qualifiers["gene"] if 'gene' in new_feature.qualifiers else new_feature.qualifiers.get("product", [f"gene {len(new_features)+1}"])
            new_features.append(new_feature)


    reference.features = new_features
    with open(f"{args.output_dir}/files/genemap.gff", "w") as f:
        GFF.write([reference], f)
