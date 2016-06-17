Introduction:

SampleSpecificDBGenerator is a program that takes RNA sequencing data analysis results,
and translates them into protein database entries that are appended to a protein database.
More information on the data analysis files can be found below. The protein database can 
either be in UniProt-XML format (by default) or protein fasta format(using the -z option).
The UniProt-XML is condensed, meaning that much of the extraneous author information is 
removed to speed up search times in Morpheus.


Example commands:

1. Create a sample-specific UniProt-XML database from a reference XML:
python samplespecificdbgenerator.py -x uniprot.xml -p ensembl.pep.all.fasta 
-g ensembl.gtf -v snpeff.vcf -b tophat.splice.bed -o sample_specific.xml

2. Sample-specific FASTA from a reference protein FASTA:
python samplespecificdbgenerator.py -p ensembl.pep.all.fasta 
-g ensembl.gtf -v snpeff.vcf -b tophat.splice.bed -o sample_specific.fasta -z

3. Sample-specific UniProt-XML from a reference protein FASTA:
python samplespecificdbgenerator.py -p ensembl.pep.all.fasta 
-g ensembl.gtf -v snpeff.vcf -b tophat.splice.bed -o sample_specific.fasta

4. Sample-specific FASTA from a reference XML:
python samplespecificdbgenerator.py -x uniprot.xml -p ensembl.pep.all.fasta 
-g ensembl.gtf -v snpeff.vcf -b tophat.splice.bed -o sample_specific.fasta -z


Author information: Anthony Cesnik, UW-Madison


Please cite:

This program was published in:
4. Cesnik, et al. "Human Proteomic Variation Revealed by Combining RNA-Seq 
Proteogenomics and Global Post-Translational Modification (G-PTM) Search 
Strategy. J. Proteome Res. 2016, 15, 800–808.

And it is based on the following papers:
1. Sheynkman, et al. "Discovery and Mass Spectrometric Analysis of Novel 
Splice-Junction Peptides Using RNA-Seq." Mol Cell Proteomics 2013, 12, 
2341-2353.
2. Sheynkman, et al. "Large-scale mass spectrometric detection of variant
peptides resulting from nonsynonymous nucleotide differences." J Proteome 
Research 2014, 13, 228-240.
3. Sheynkman, et al. "Using Galaxy-P to leverage RNA-Seq for the discovery 
of novel protein variations." BMC Genomics 2014, 15, 9.


System requirements:
- 8 GB of RAM is recommended
- python v2.7.10
See https://www.python.org/downloads/ for installation instructions.
This includes the “pip” package manager.
- Biopython python package
Install using the command: pip install biopython
Or see http://biopython.org/wiki/Download for installation instructions.
- Lxml python package
Install using the command: pip install lxml 
Or see http://lxml.de for installation instructions.
- If you encounter errors installing either package, we recommend
trying an alternate package manager, such as Canopy, which can be found
here: https://www.enthought.com/products/canopy/.

Version updates:
v0.0.2 November 26, 2015 Initial commit
v0.0.3 November 30, 2015 Updated usage information. Allows minimum length 
cutoff to filter both SAV and NSJ peptide entries.


Usage: samplespecificdbgenerator.py [options]

Options:
  -h, --help            show this help message and exit
  -x REFERENCE_XML, --reference_xml=REFERENCE_XML
                        Reference protein UniProt-XML file. Sequence variant
                        peptide entries are appended to this database to
                        generate the ouptut UniProt-XML protein database.
  -p PROTEIN_FASTA, --protein_fasta=PROTEIN_FASTA
                        Reference protein FASTA file. Used to generate SAV
                        peptide entries. If no UniProt-XML is specified, SAV
                        and NSJ entries will be appended to this database to
                        generate an output database. By default, this output
                        will be a UniProt-XML protein database without PTM
                        annotations. If --output-fasta is selected, the output
                        will be a protein FASTA.
  -g GENE_MODEL, --gene_model=GENE_MODEL
                        GTF gene model file. Used to annotate NSJ peptide
                        entries.
  -v SNPEFF_VCF, --snpeff_vcf=SNPEFF_VCF
                        SnpEff VCF file with HGVS annotations (else read from
                        stdin).
  -b SPLICE_BED, --splice_bed=SPLICE_BED
                        BED file (tophat junctions.bed) with sequence column
                        added.
  -o OUTPUT, --output=OUTPUT
                        Output file path. Outputs UniProt-XML format unless
                        --output-fasta is selected.
  -z, --output_fasta    Output a FASTA-format database. Place path for output
                        file after the --output flag.
  -l LEADING_AA_NUM, --leading_aa_num=LEADING_AA_NUM
                        Leading number of AAs to output for SAV peptides.
                        Default: 33.
  -t TRAILING_AA_NUM, --trailing_aa_num=TRAILING_AA_NUM
                        Trailing number of AAs to output for SAV peptides.
                        Default: 33.
  -D NSJ_DEPTH_CUTOFF, --nsj_depth_cutoff=NSJ_DEPTH_CUTOFF
                        Keep only NSJs found with above this depth (BED score
                        field). Default: 0.
  -E SNV_DEPTH_CUTOFF, --snv_depth_cutoff=SNV_DEPTH_CUTOFF
                        Keep only SNVs found with above this depth (DP=#
                        field). Default: 0.
  -M MINIMUM_LENGTH, --minimum_length=MINIMUM_LENGTH
                        Keep only sequence variant peptides with greater than
                        or equal to this length. Default: 0.
  -Q BED_SCORE_NAME, --bed_score_name=BED_SCORE_NAME
                        Include in the NSJ ID line score_name:score. Default:
                        "depth."
  -R REFERENCE, --reference=REFERENCE
                        Genome Reference Name for NSJ ID location.
                        Automatically pulled from genome_build header in GTF
                        if present.
                        