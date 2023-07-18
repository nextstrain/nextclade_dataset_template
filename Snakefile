build_dir = config.get('build_dir', 'results')
genes = config.get('genes', [])

if "clades" in config:
    print(f"using clades-tsv '{config['clades']}' to assign clades")
    ruleorder: clades_from_tsv>clades_from_metadata
else:
    print(f"using metadata to assign clades")
    ruleorder: clades_from_metadata>clades_from_tsv

rule align:
    input:
        sequences = config["sequences"],
        ref = config["reference"],
        annotation = config["annotation"]
    output:
        alignment = build_dir + "/alignment.fasta",
        translations = expand(build_dir + "/translation_{gene}.fasta", gene=genes),
        insertions = build_dir + "/insertions.csv"
    params:
        translations = lambda w:build_dir + "/translation_{gene}.fasta"
    threads: 4
    shell:
        """
        ./nextalign run -j {threads} -r {input.ref} -m {input.annotation} --output-fasta {output.alignment} --output-translations {params.translations} \
                            --output-insertions {output.insertions} --include-reference --gap-alignment-side right \
                        --penalty-gap-open 12 --penalty-gap-open-out-of-frame 14 --penalty-gap-open-in-frame 10 \
                        {input.sequences}
        """

rule tree:
    input:
        aln = build_dir + "/alignment.fasta"
    output:
        tree = build_dir + "/tree_raw.nwk"
    shell:
        """
        augur tree --alignment {input.aln} --output {output.tree}
        """

rule refine:
    input:
        tree = build_dir + "/tree_raw.nwk",
        metadata = config["metadata"],
        alignment = build_dir + "/alignment.fasta"
    output:
        tree = build_dir + "/tree.nwk",
        node_data = build_dir + "/branch_lengths.json"
    params:
        options = config.get("refine_options", "")
    shell:
        """
        augur refine --tree {input.tree} --metadata {input.metadata} \
                     --alignment {input.alignment} {params.options} \
                     --output-tree {output.tree} --output-node-data {output.node_data}
        """


rule ancestral:
    input:
        tree = build_dir + "/tree.nwk",
        aln = build_dir + "/alignment.fasta"
    output:
        node_data=build_dir + "/nt_muts.json"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.aln} --output-node-data {output.node_data}
        """

rule translate:
    input:
        tree = build_dir + "/tree.nwk",
        ref = config["reference"],
        aln = build_dir + "/alignment.fasta"
    output:
        node_data = build_dir + "/aa_muts.json"
    params:
        translations = expand(build_dir + "/translation_{gene}.fasta", gene=genes),
        genes = genes
    shell:
        """
        python3 scripts/explicit_translation.py --tree {input.tree} --reference {input.ref}\
                    --translations {params.translations} --genes {params.genes}\
                    --output {output.node_data}
        """

rule clades_from_metadata:
    input:
        tree = build_dir + "/tree.nwk",
        metadata = config["metadata"]
    output:
        build_dir + "/clades.json"
    params:
        clade_key = config.get("clade_key", "clade"),
    shell:
        """
        python3 scripts/assign_clades.py --tree {input.tree} --metadata {input.metadata} \
                                         --clade-key {params.clade_key} --output {output}
        """

rule clades_from_tsv:
    input:
        tree = build_dir + "/tree.nwk",
        clades = config.get("clades",''),
        muts = [build_dir + "/nt_muts.json", build_dir + "/aa_muts.json"]
    output:
        build_dir + "/clades.json"
    shell:
        """
        augur clades --tree {input.tree} --mutations {input.muts} --clades {input.clades} \
                     --output {output}
        """


# make sure all differences between the alignment reference and the root are attached as mutations to the root
rule attach_root_mutations:
    input:
        aa_muts=build_dir + "/aa_muts.json",
        nuc_muts=build_dir + "/nt_muts.json",
        tree = build_dir + "/tree.nwk"
    output:
        aa_muts=build_dir + "/aa_muts_adapted.json",
        nuc_muts=build_dir + "/nt_muts_adapted.json"
    params:
        genes = genes,
        translations= expand(build_dir + "/translation_{gene}.fasta",gene=genes),
        reference = config['reference_name']
    shell:
        """
        python3 scripts/attach_root_mutations.py \
            --tree {input.tree} \
            --translations {params.translations:q} \
            --reference {params.reference} \
            --genes {params.genes} \
            --aa-mutations {input.aa_muts} \
            --nuc-mutations {input.nuc_muts} \
            --output-aa-mutations {output.aa_muts} \
            --output-nuc-mutations {output.nuc_muts}
        """

rule export:
    input:
        tree = build_dir + "/tree.nwk",
        node_data = [build_dir + "/branch_lengths.json", build_dir + "/clades.json",
                     build_dir + "/nt_muts_adapted.json", build_dir + "/aa_muts_adapted.json"],
        auspice_config = config["auspice_config"]
    output:
        build_dir + "/tree.json"
    shell:
        """
        augur export v2 --tree {input.tree} --node-data {input.node_data} \
                        --auspice-config {input.auspice_config} --output {output} --minify-json
        """


rule assemble:
    input:
        tree = build_dir + "/tree.json",
        reference = config["reference"],
        annotation = config["annotation"],
        primers = config.get("primers", "dummy_files/primers.csv"),
        virus_properties = config.get("virus_properties", "dummy_files/virus_properties.json"),
        qc = config.get("qc_json", "dummy_files/qc.json"),
        tag = config.get("tag_json", "dummy_files/tag.json")
    output:
        config['dataset_dir'] + "/tree.json",
    params:
        dataset_dir = config['dataset_dir'],
        example_data = config.get("example_data", "")
    shell:
        """
        cp {input} {params.example_data} {params.dataset_dir}
        """
