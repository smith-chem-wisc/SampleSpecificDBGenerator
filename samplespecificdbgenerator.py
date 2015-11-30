# main program: samplespecificdbgenerator.py

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:41:36 PM$"

import sys
import os.path
import optparse
import refparse
import novelsplices
import variantcalls
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'
    
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
      #I/O
    parser.add_option( '-x', '--reference_xml', dest='reference_xml', help='Reference protein UniProt-XML file. Sequence variant peptide entries are appended to this database to generate the ouptut UniProt-XML protein database.' )
    parser.add_option( '-p', '--protein_fasta', dest='protein_fasta', help='Reference protein FASTA file. Used to generate SAV peptide entries. If no UniProt-XML is specified, SAV and NSJ entries will be appended to this database to generate an output database. By default, this output will be a UniProt-XML protein database without PTM annotations. If --output-fasta is selected, the output will be a protein FASTA.')
    parser.add_option( '-g', '--gene_model', dest='gene_model', default=None, help='GTF gene model file. Used to annotate NSJ peptide entries.')
    parser.add_option( '-v', '--snpeff_vcf', dest='snpeff_vcf', help='SnpEff VCF file with HGVS annotations (else read from stdin).' )
    parser.add_option( '-b', '--splice_bed', dest='splice_bed', help='BED file (tophat junctions.bed) with sequence column added.' )
    parser.add_option( '-o', '--output', dest='output', help='Output file path. Outputs UniProt-XML format unless --output-fasta is selected.' )
    parser.add_option( '-z', '--output_fasta', dest='output_fasta', action='store_true', default=False, help='Output a FASTA-format database. Place path for output file after the --output flag.')
      #Peptide sequence construction
    parser.add_option( '-l', '--leading_aa_num', dest='leading_aa_num', type='int', default=33, help='Leading number of AAs to output for SAV peptides. Default: 33.' )
    parser.add_option( '-t', '--trailing_aa_num', dest='trailing_aa_num', type='int', default=33, help='Trailing number of AAs to output for SAV peptides. Default: 33.' )
      #Filtering parameters
    parser.add_option( '-D', '--nsj_depth_cutoff', dest='nsj_depth_cutoff', type='int', default=0, help='Keep only NSJs found with above this depth (BED score field). Default: 0.' )
    parser.add_option( '-E', '--snv_depth_cutoff', dest='snv_depth_cutoff', type='int', default=0, help='Keep only SNVs found with above this depth (DP=# field). Default: 0.' )
    parser.add_option( '-M', '--minimum_length', dest='minimum_length', type='int', default=0, help='Keep only sequence variant peptides with greater than or equal to this length. Default: 0.' )
      #Simple entry
    parser.add_option( '-Q', '--bed_score_name', dest='bed_score_name', default="depth", help='Include in the NSJ ID line score_name:score. Default: "depth."'  )
    parser.add_option( '-R', '--reference', dest='reference', default="None", help='Genome Reference Name for NSJ ID location. Automatically pulled from genome_build header in GTF if present.'  )
    (options, args) = parser.parse_args()
    
    ##INPUTS##
    #Protein FASTA
    try:
        protein_fasta = os.path.abspath(options.protein_fasta)
        protein_fasta = open(protein_fasta, 'r')
        protein_fasta = refparse.read_protein_fasta(protein_fasta)
    except Exception, e:
        print >> sys.stderr, "failed: %s" % e
        exit(2)
    
    #Reference XML/FASTA
    if options.reference_xml != None:
        try:
            refXml = os.path.abspath(options.reference_xml)
            refXml = open(refXml, 'r')
            p = et.XMLParser(remove_blank_text=True) #required for pretty additions
            db = et.parse(refXml, p)
            root = db.getroot()
            for entry in root:
                refparse.condense_xml_entry(entry)
        except Exception, e:
            print >> sys.stderr, "Parsing and/or condensing reference xml failed: %s" % e
            exit(2)
    elif options.protein_fasta != None:
        try:
            refFasta = os.path.abspath(options.protein_fasta)
            refFasta = open(refFasta, 'r')
            root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
            db = et.ElementTree(root)
            refparse.read_fasta_to_xml(root, refFasta)
        except Exception, e:
            print >> sys.stderr, "Parsing reference fasta failed: %s" %e
            exit(2)
    else:
        print >> sys.stderr, "failed: no reference database specified"

    ##OUTPUT##
    outFile = None
    if options.output == None: outFile = sys.stdout
    else:
        try:
            outFile = os.path.abspath(options.output)
            outFile = open(outFile, 'w')
        except Exception, e:
            print >> sys.stderr, "Opening outfile failed: %s" % e
            exit(3)
            
    #Process gene model
    try:
        geneModelFile = os.path.abspath(options.gene_model)
        linect = sum(1 for line in open(geneModelFile))
        geneModelFile = open(geneModelFile, 'r')
        geneModel = refparse.GeneModel()
        for i, gtf_line in enumerate(geneModelFile):
            if i % 20000 == 0: print "gene_model line " + str(i) + " of " + str(linect)
            if gtf_line.startswith('#') and  gtf_line.find('genome_build') >= 0: options.reference = gtf_line.split()[1]
            elif gtf_line.startswith('#'): continue
            else: geneModel.new_entry(gtf_line)
    except Exception, e:
        print >> sys.stderr, "Parsing gene model failed: %s" % e
        exit(2)
    
    #Process VCF
    try:
        snpeff_vcf = os.path.abspath(options.snpeff_vcf)
        linect = sum(1 for line in open(snpeff_vcf))
        snpeff_vcf = open(snpeff_vcf, 'r')
        for i, vcf_line in enumerate(snpeff_vcf):
            if i % 1000 == 0: print "snpeff_vcf line " + str(i) + " of " + str(linect)
            variantcalls.parse_vcf_line(root, vcf_line, protein_fasta, options.snv_depth_cutoff, options.minimum_length, options.leading_aa_num, options.trailing_aa_num)
        snpeff_vcf.close()
    except Exception, e:
        print >> sys.stderr, "VCF processing failed: %s" % e
        exit(1)
    
    #Process Splice BED
    try:
        splice_bed = os.path.abspath(options.splice_bed)
        linect = sum(1 for line in open(splice_bed))
        splice_bed = open(splice_bed, 'r')
        for i, bed_line in enumerate(splice_bed):
            if i % 2000 == 0: print "splice_bed line " + str(i) + " of " + str(linect)
            if bed_line.startswith('track'): continue
            novelsplices.translate_bed_line(root, geneModel, bed_line, options.nsj_depth_cutoff, options.minimum_length, options.reference, options.bed_score_name)
        splice_bed.close()
    except Exception, e:
        print >> sys.stderr, "Splice BED processing failed: %s" % e
        exit(2)
    
    #Write the sample specific database to outfile
    if not options.output_fasta: db.write(outFile, pretty_print=True)
    else:
        entryct = len(root)
        if outFile != None:
            for i, entry in enumerate(root):
                if i % 1000 == 0: print "writing entry " + str(i) + " of " + str(entryct)
                entry = refparse.xml_to_fasta(entry)
                if entry == None: continue
                else: outFile.write(entry[0] + '\n' + entry[1] + '\n')
            outFile.close()
        else: print >> sys.stderr, "Writing to fasta failed: no out-file found"

if __name__ == "__main__" : __main__()

