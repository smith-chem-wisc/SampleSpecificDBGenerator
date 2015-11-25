# test module: test_variantcalls.py
# main program: samplespecificdbgenerator.py

import unittest
import variantcalls
import refparse
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

class  Test_enter_seq_var(unittest.TestCase):
    def setUp(self):
        self.root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
        self.db = et.ElementTree(self.root)

    def test_enter_seqvar(self):
        variantcalls.enter_seqvar(self.root, ">accession", "peptide:seqtype", "chromosome", "full name", "loci", "score", "geneId", "transcriptId", "seq\nuence")
        self.assertNotEqual(self.root[0].find(UP+'accession'), None)
        self.assertTrue(self.root[0].find(UP+'accession').text not in [None, ''])
        self.assertNotEqual(self.root[0].find(UP+'name'), None)
        self.assertTrue(self.root[0].find(UP+'name').text not in [None, ''])
        self.assertNotEqual(self.root[0].find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName'), None)
        self.assertTrue(self.root[0].find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName').text not in [None, ''])
        self.assertNotEqual(self.root[0].find(UP+'organism'), None)
        self.assertNotEqual(self.root[0].find(UP+'proteinExistence'), None)
        self.assertNotEqual(self.root[0].find(UP+'proteinExistence').find(UP+'depth'), None)
        self.assertNotEqual(self.root[0].find(UP+'sequence'), None)
        self.assertTrue(self.root[0].find(UP+'sequence').text not in [None, ''])
        self.assertTrue(self.root[0].find(UP+'sequence').text.find('\n') < 0)

class Test_parse_aa_change(unittest.TestCase):
    def test_parse_aa_change_standard_format(self):
        aa_change = "G528R"
        aa_pos, ref_aa, alt_aa = variantcalls.parse_aa_change(aa_change)
        self.assertTrue(aa_pos == 528)
        self.assertTrue(ref_aa == 'G')
        self.assertTrue(alt_aa == 'R')
    def test_parse_aa_change_alternate_format(self):
        aa_change = "p.Gly528Arg/c.1582G>C"
        aa_pos, ref_aa, alt_aa = variantcalls.parse_aa_change(aa_change)
        self.assertTrue(aa_pos == 528)
        self.assertTrue(ref_aa == 'G')
        self.assertTrue(alt_aa == 'R')
        
class Test_transcript_based_entry(unittest.TestCase):
    def setUp(self):
        self.root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
        self.db = et.ElementTree(self.root)
        fastaLines = ""
        fastaLines += ">ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding"
        fastaLines += "\nMPFLELDTNLPANRVPAGLEKRLCAAAASILGKPADRVNVTVRPGLAMALSGSTEPCAQL"
        fastaLines += "\nSISSIGVVGTAEDNRSHSAHFFEFLTKELALGQDRILIRFFPLESWQIGKIGTVMTFL"
        self.header = ">ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding"        
        self.sequence = "MPFLELDTNLPANRVPAGLEKRLCAAAASILGKPADRVNVTVRPGLAMALSGSTEPCAQLSISSIGVVGTAEDNRSHSAHFFEFLTKELALGQDRILIRFFPLESWQIGKIGTVMTFL"
        self.protein_fasta = open('fasta_entry.fasta', 'w')
        self.protein_fasta.write(fastaLines)
        self.protein_fasta.close()
        self.protein_fasta = open('fasta_entry.fasta', 'r')
        self.protein_fasta = refparse.read_protein_fasta(self.protein_fasta)
    
    def test_transcript_based_entry(self):
        variantcalls.transcript_based_entry(self.root, '', "ENST00000398344", self.protein_fasta, "22", "24313554:24316773", "codon_change", "A12B", "B", 12, 0, 33, 33)
        self.assertTrue(self.root[0].find('.//'+UP+'fullName').text.find('sav:A12B') >= 0)
        self.assertTrue(self.root[0].find(UP+'sequence').text.find('B') >= 0)
        self.assertTrue(self.root[0].find(UP+'sequence').text[11] == 'B') #12 is 1-indexed; 11 is 0-indexed

if __name__ == '__main__':
    unittest.main()

