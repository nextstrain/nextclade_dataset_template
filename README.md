## Workflow to boot-strap nextclade datasets

Nextclade requires datasets for pathogen specific information and rich datasets allow sophisticated analysis.
But many of the basic functionalities like alignment and mutations calling only require a reference sequence and its annotation.


### Generate a nextclade conforming annotation

This step requires some user interaction. You need to specify the accession number of a well annotated sequence, ideally from RefSeq, and run
```
python scripts/generate_from_genbank.py --reference <accession> --output-dir <output-directory>
```
The script will fetch the sequence and annotation from NCBI and will ask you to rename CDS or proteins to sensible names. This interactive procedure is only sensible for small viral genomes, for larger genomes like MPXV you will want to use a different method.
The scripts outputs a `reference.fasta` and an `annotation.gff` in the specified directory.
These two files are required to align and translate using `nextclade` and `nextalign`.
We suggest putting these files into the profile associated with the dataset.

For HIV-1, this could look like this:
```
python scripts/generate_from_genbank.py --reference NC_001802 --output-dir profiles/HIV-1
```


### Input files required for clade annotation and QC
To allow nextclade to annotate clades and determine reversion mutations, private mutations, and other diagnostic features, you need a reference tree and a qc config.
The tree can be generated from a set of high quality genomes and some associated metadata.








