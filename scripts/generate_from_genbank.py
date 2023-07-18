#!/bin/env python3
import shutil
from Bio import SeqIO
from collections import defaultdict


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

def get_gff(accession):
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={accession}"
    import urllib
    return [x.decode() for x in urllib.request.urlopen(url).readlines()]


if __name__=="__main__":
    args = parse_args()
    reference = get_reference_sequence(args.reference)
    gff = get_gff(args.reference)

    # write the reference sequence to the output directory
    SeqIO.write(reference, f"{args.output_dir}/reference.fasta", "fasta")

    # collect all cds and alternative annotations of protein sequences
    all_cds = defaultdict(lambda: defaultdict(list))
    for line in gff:
        if line[0]=='#':
            continue
        entries = line.strip().split('\t')
        if len(entries)<9:
            continue # invalid line

        # parse attributes and feature type and ID
        attributes = {x.split('=')[0]:x.split('=')[1] for x in entries[-1].split(';')}
        feature_type = entries[2]
        # IDs look like this: ID=id-NP_057850.1:133..363 where the part after to colon is the range in the translated sequence (exists in cases of mature protein annotations)
        feature_id = attributes['ID'].split(':')[0].split('-')[-1]
        if feature_type=='CDS':
            all_cds[feature_id]['CDS'].append([entries[:-1], attributes])
        elif feature_type in ['mature_protein_region_of_CDS', 'mat_peptide', 'mat_protein']:
            all_cds[feature_id]['mature_protein'].append([entries[:-1], attributes])


    streamlined_cds = {}
    names_by_id = {}
    # loop through all CDS and ask the user to pick one of the annotations
    # allow renaming of the CDS to user friendly names
    for cds_id, cds_sets in all_cds.items():
        if len(cds_sets)>1:
            print(f"\nCDS {cds_id} is annotated in multiple ways:")
            for i, (cds_set, segments) in enumerate(cds_sets.items()):
                print(f"\t{cds_set} with a total of {len(segments)} items [{i+1}]")
            choice = int(input("Please pick number in brackets to choose one (0 for omission): "))
            if choice:
                cds_set = list(cds_sets.keys())[choice-1]
                segments = cds_sets[cds_set]
            else:
                continue
        else:
            cds_set = list(cds_sets.keys())[0]
            segments = cds_sets[cds_set]

        for segment in segments:
            segment_id = segment[1]["ID"]
            # only needed for one segment per CDS, skip if already in streamlined_cds
            if segment_id not in streamlined_cds:
                streamlined_cds[segment_id] = []

                print(f"Attributes of the segment with ID='{segment_id}' are:")
                for k,v in segment[1].items():
                    print(f'\t\t{k}:\t{v}')
                new_name = input("Enter desired name (leave empty to drop): ")
                if new_name:
                    names_by_id[segment_id] = new_name
                else:
                    continue
            # if renamed and selected, add the segment to the list of segments
            if segment_id in names_by_id:
                new_entries, new_attributes = list(segment[0]), dict(segment[1])
                new_entries[2]='CDS'
                new_attributes['Name']=names_by_id[segment_id]
                if "Parent" in new_attributes: new_attributes.pop("Parent")
                streamlined_cds[segment_id].append([new_entries, new_attributes])

    # write the gff file as a simple text file line by line
    with open(f"{args.output_dir}/genemap.gff", "w") as f:
        # write the header and the region line
        for line in gff:
            entries = line.split('\t')
            if entries[0][0]=='#' or (len(entries)>8 and entries[2]=='region'):
                f.write(line)

        # write the CDS lines
        for cds in streamlined_cds:
            for segment in streamlined_cds[cds]:
                attributes = ';'.join([f"{k}={v}" for k,v in sorted(segment[1].items(), key=lambda x:(x[0]!='Name', len(x[1])))])
                f.write('\t'.join(segment[0])+'\t' + attributes + '\n')

