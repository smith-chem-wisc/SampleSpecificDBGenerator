# module name: variantcalls
# main program: samplespecificdbgenerator
#
# Amino_Acid_Change notations
# G528R
# p.Gly528Arg/c.1582G>C
#
# This is the current format of the EFF entry:
# EFF=missense(MODERATE|MISSENSE|Ggg/Cgg|G528R|802|SCNN1D|protein_coding|CODING|ENST00000379116|12|1);OICR=(ENST00000379116|1808)
# If this becomes variable, will need to dynamically pattern this on the defintion in the vcf header:
# #INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank | Genotype_Number [ | ERRORS | WARNINGS ] )' ">
#
# Header:
# >ENSP00000317992 pep:sav chromosome:GRCh37:1:879584:894670:-1 gene:ENSG00000188976 transcript:ENST00000327044 gene_biotype:protein_coding transcript_biotype:protein_coding snv_location:1:888659 codon_change:Atc/Gtc sav:I300V

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:55:09 PM$"

import sys
import re
from lxml import etree as et
import refparse

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def enter_seqvar(root, acc, seqtype, chromosome, addedInfo, loci, score, geneId, transcriptId, seq):
    entry = et.SubElement(root, UP+'entry', dataset="Ensembl")

    accession = et.SubElement(entry, UP+'accession')
    accession.text = acc.strip('>')

    #Set the names
    name = et.SubElement(entry, UP+'name')
    name.text = seqtype

    fullName = et.SubElement(et.SubElement(et.SubElement(entry, UP+'protein'), UP+'recommendedName'), UP+'fullName')
    fullName.text = addedInfo.replace('\n','').replace('\r','')

    gene = et.SubElement(entry, UP+'gene')
    geneName1 = et.SubElement(gene, UP+'name', type="coords")
    if loci: geneName1.text = chromosome + " " + loci
    else: geneName1.text = chromosome
    geneName2 = et.SubElement(gene, UP+'name', type="primary")
    geneName2.text = geneId
    transcript = et.SubElement(geneName2, UP+'transcript', type="primary")
    transcript.text = transcriptId

    organism = et.SubElement(entry, UP+'organism')
    organismName1 = et.SubElement(organism, UP+'name', type="scientific")
    organismName2 = et.SubElement(organism, UP+'name', type="common")
    organismName1.text = "Homo sapiens"
    organismName2.text = "Human"

    proteinExist = et.SubElement(entry, UP+'proteinExistence', type="evidence at transcript level")
    if score: et.SubElement(proteinExist, UP+'depth', reads=str(score))

    sequence = et.SubElement(entry, UP+'sequence', version="1", fragment="single")
    sequence.text = seq.replace('\n','').replace('\r','')

def parse_aa_change(aa_change):
    aa_abbrev_dict = refparse.aa_abbrev_dict()
    aa_change_regex = '([A-Z])(\d+)([A-Z])' # G528R
    aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))' # p.Gly528Arg/c.1582G>C
    aa_pos = None # 1-based position
    ref_aa, alt_aa = '_', '_'
    m = re.match(aa_change_regex, aa_change) # parse aa_change, and get AA change position and alternate Animo Acid
    if m:
        aa_pos = int(m.groups()[1])
        ref_aa = m.groups()[0]
        alt_aa = m.groups()[2]
    else:
        m = re.match(aa_hgvs_regex, aa_change)
        if m:
            aa_pos = int(m.groups()[1])
            ref_aa = aa_abbrev_dict[m.groups()[0]]
            alt_aa = aa_abbrev_dict[m.groups()[2]]
    return aa_pos, ref_aa, alt_aa
            
def transcript_based_entry(root, line, transcript, protein_fasta, chrom, pos, codon_change, sav, alt_aa, aa_pos, minPepLength, leading_aas, trailing_aas):
    (pep_id, pep_seq) = refparse.get_protein_fasta_seq(transcript, protein_fasta) # get AA sequence using transcript ID in VCF line
    aa_offset = (aa_pos - 1)
    start_pos, end_pos = 0, 0 
    if not pep_seq: 
        print "error finding transcript " + transcript + " in " + line
        return
        
    start_pos = max(aa_offset - leading_aas, 0) if leading_aas else 0
    end_pos = min(aa_offset + trailing_aas + 1, len(pep_seq)) if trailing_aas else len(pep_seq)
    alt_seq = pep_seq[start_pos:aa_offset] + alt_aa + pep_seq[aa_offset+1:end_pos]
    if len(alt_seq) < minPepLength: return

    fields = pep_id.split()
    fields.append("snv_location:%s:%s codon_change:%s sav:%s" % (chrom, pos, codon_change, sav))
    enter_seqvar(root, fields[0], 'pep:sav', fields[2], ' '.join(fields[5:]), '', '', fields[3][5:], fields[4][11:], alt_seq)
    
def parse_vcf_line(root, line, protein_fasta, snvDepthCut, minPepLength, leading_aas, trailing_aas):
    vcf_header = [] # Save VCF file header, not currently used
    if line.startswith('##'): vcf_header.append(line) ##SnpEffVersion #May need to check SnpEff version in the header, the EFF info changed between versions 2 and 3
    elif line.startswith('#CHROM'): return
    else:
        fields = line.split('\t')
        (chrom, pos, id, ref, alts, qual, filter, info) = fields[0:8]
        qual = float(qual)
        for info_item in info.split(';'):
            try:
                if info_item.find('=') < 0: return
                (key, val) = info_item.split('=', 1)
                if key == 'DP' and int(val) <= snvDepthCut: return # filter out the entries that have poor coverage AC150310
                if key == 'EFF':
                    for effect in val.split(','):
                        (eff, effs) = effect.rstrip(')').split('(')                                
                        if eff not in ['NON_SYNONYMOUS_CODING','MISSENSE','missense_variant']: continue #updated for snpeff.4.0 with the inclusion of MISSENSE
                        (impact, functional_class, codon_change, aa_change, aa_len, gene_name, biotype, coding, transcript, exon) = effs.split('|')[0:10]
                        if transcript:
                            aa_pos, ref_aa, alt_aa = parse_aa_change(aa_change)
                            if not aa_pos: continue
                            sav = "%s%d%s" % (ref_aa, aa_pos, alt_aa)
                            transcript_based_entry(root, line, transcript, protein_fasta, chrom, pos, codon_change, sav, alt_aa, aa_pos, minPepLength, leading_aas, trailing_aas)
            except Exception, e:
                print >> sys.stderr, "failed: %s" % e
