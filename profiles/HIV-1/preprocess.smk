subtypes = ["B", "C", "01_AE", "A1", "A1D", "01B", "D",
            "A1C", "02_AG", "BF1", "A1CD", "A6", "CD",
            "BC", "01BC", "G", "0107", "F1", "07BC", "08_BC"]


rule parse:
    input:
        sequences = "data/HIV-1/2023-07-18_LANL_download.fasta.xz",
        reference = "profiles/HIV-1/reference.fasta"
    output:
        sequences = "data/HIV-1/2023-07-18_parsed.fasta.xz",
        metadata= "data/HIV-1/2023-07-18_metadata.tsv"
    params:
        separator = '.',
        non_acgt_count = 30
    run:
        from Bio import SeqIO, SeqRecord
        import lzma
        import pandas as pd
        metadata = [{'strain':'NC_001802','subtype':'B', 'country':'FR',
                     'date':"1983-XX-XX", 'accession':'NC_001802', 'name':"HXB2_reference"}]
        fields = {'strain': -1, 'subtype': 0, 'country': 1, 'date': 2, 'accession': -1,  'name': [3,-1]}
        fh_out = lzma.open(output.sequences, "wt")
        SeqIO.write(SeqIO.read(input.reference, 'fasta'), fh_out, "fasta")
        with lzma.open(input.sequences, "rt") as fh_in:
            for record in SeqIO.parse(fh_in, "fasta"):
                entries = record.id.split('.')
                if len(entries)<5:
                    print(entries)
                    continue
                datum = {k: entries[v] if type(v)==int else '.'.join(entries[v[0]:v[1]]) for k, v in fields.items()}
                if datum['subtype'] not in subtypes:
                    continue
                try:
                    datum['date'] = f"{int(datum['date'])}-XX-XX"
                except:
                    datum['date'] = None
                record.id = datum['strain']
                record.seq = record.seq.upper()
                non_acgt = len(record.seq) - sum([record.seq.count(x) for x in 'ACGT-'])
                if non_acgt<params.non_acgt_count:
                    SeqIO.write(record, fh_out, "fasta")
                    metadata.append(datum)
                else:
                    print(record.id, 'excluded', non_acgt)

            pd.DataFrame(metadata).to_csv(output.metadata, sep='\t', index=False)
        fh_out.close()

rule filter:
    input:
        sequences = "data/HIV-1/2023-07-18_parsed.fasta.xz",
        reference = "profiles/HIV-1/reference.fasta",
        metadata= "data/HIV-1/2023-07-18_metadata.tsv"
    output:
        sequences = "data/HIV-1/sequences.fasta",
        metadata= "data/HIV-1/metadata.tsv"
    params:
        group_by = 'subtype year',
        max_number_sequences = 500
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.max_number_sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """



rule sliced_alignment_for_tree:
    input:
        build_dir + "/alignment.fasta"
    output:
        config.get("tree_alignment", 'tmp')
    run:
        from Bio import SeqIO
        with open(output[0], "w") as f:
            for record in SeqIO.parse(input[0], "fasta"):
                record.seq = record.seq[1798:3415]
                SeqIO.write(record, f, "fasta")
