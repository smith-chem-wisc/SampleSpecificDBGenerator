# test module: test_refparse.py
# main program: samplespecificdbgenerator.py

import unittest
import refparse
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

class  Test_condense_xml_entry(unittest.TestCase):
    def setUp(self):
        self.root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
        self.db = et.ElementTree(self.root)
    
    def test_condense_xml_entry(self):
        #read in rawuniprotentry.xml
        for entry in self.root:
            refparse.condense_xml_entry(entry)
            self.assertNotEqual(self.entry.find(UP+'accession'), None)
            self.assertTrue(self.entry.find(UP+'accession').text not in [None, ''])
            self.assertNotEqual(self.entry.find(UP+'name'), None)
            self.assertTrue(self.entry.find(UP+'name').text not in [None, ''])
            self.assertNotEqual(self.entry.find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName'), None)
            self.assertTrue(self.entry.find(UP+'protein').find(UP+'recommendedName').find(UP+'fullName').text not in [None, ''])
            self.assertNotEqual(self.entry.find(UP+'organism'), None)
            self.assertNotEqual(self.entry.find(UP+'proteinExistence'), None)
            self.assertNotEqual(self.entry.find(UP+'proteinExistence').find(UP+'depth'), None)
            self.assertNotEqual(self.entry.find(UP+'sequence'), None)
            self.assertTrue(self.entry.find(UP+'sequence').text not in [None, ''])
        
class Test_fasta_mainuplation(unittest.TestCase):
    def setUp(self):
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
        
    def test_read_protein_fasta(self):
        self.assertTrue(self.sequence in self.protein_fasta[1])
        self.assertTrue(self.protein_fasta[0][0].find("ENST00000398344") >= 0)
    
    def test_get_protein_fasta_seq(self):
        self.assertTrue(self.header, self.sequence == refparse.get_protein_fasta_seq("ENST00000398344", self.protein_fasta))
    
class Test_parse_gtf_line(unittest.TestCase):
    def test_parse_gtf_line(self):
        line1="1\tprocessed_transcript\texon\t12613\t12721\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\"; exon_number \"2\"; gene_name \"DDX11L1\"; gene_biotype \"pseudogene\"; transcript_name \"DDX11L1-002\"; exon_id \"ENSE00003582793\";"
        line2="MT\tprotein_coding\texon\t14149\t14673\t.\t-\t.\tgene_id \"ENSG00000198695\"; transcript_id \"ENST00000361681\"; exon_number \"1\"; gene_name \"MT-ND6\"; gene_biotype \"protein_coding\"; transcript_name \"MT-ND6-201\"; exon_id \"ENSE00001434974\";"
        line3="MT\tprotein_coding\texon\t14747\t15887\t.\t+\t.\tgene_id \"ENSG00000198727\"; transcript_id \"ENST00000361789\"; exon_number \"1\"; gene_name \"MT-CYB\"; gene_biotype \"protein_coding\"; transcript_name \"MT-CYB-201\"; exon_id \"ENSE00001436074\";"
        line4="MT\tprotein_coding\texon\t18000\t18002\t.\t+\t.\tgene_id \"ENSG00000198727\"; transcript_id \"ENST00000361789\"; exon_number \"1\"; gene_name \"MT-CYBz\"; gene_biotype \"protein_coding\"; transcript_name \"MT-CYB-201z\"; exon_id \"ENSE00001436074\";"
        geneModel = refparse.GeneModel()
        geneModel.new_entry(line1)
        geneModel.new_entry(line2)
        geneModel.new_entry(line3)
        geneModel.new_entry(line4)
        
        #Initializations
        self.assertTrue(('1','+') in geneModel.chromosomes) #Initialization of gene model
        self.assertTrue(('MT','+') in geneModel.chromosomes) #Initialization of second chromosome
        self.assertTrue(('MT','-') in geneModel.chromosomes) #Initialization of second strand
        
        #Boundaries of entry
        self.assertEqual(geneModel.identify_range('1', '+', 12613, 12800), ("DDX11L1","pseudogene")) #Test start-in-bounds left edge of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12700, 12800), ("DDX11L1","pseudogene")) #Test start-in-bounds central edge of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12721, 12800), ("DDX11L1","pseudogene")) #Test start-in-bounds right edge of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12612, 12800), ("DDX11L1","pseudogene")) #Test start-in-bounds spanning edges of gene model entry
        self.assertNotEqual(geneModel.identify_range('1', '+', 12722, 12800), ("DDX11L1","pseudogene")) #Test start-in-bounds not in bounds of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12612, 12613), ("DDX11L1","pseudogene")) #Test end-in-bounds left edges of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12612, 12700), ("DDX11L1","pseudogene")) #Test end-in-bounds central edges of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12612, 12721), ("DDX11L1","pseudogene")) #Test end-in-bounds right edges of gene model entry
        self.assertEqual(geneModel.identify_range('1', '+', 12612, 12722), ("DDX11L1","pseudogene")) #Test end-in-bounds spanning edges of gene model entry
        self.assertNotEqual(geneModel.identify_range('1', '+', 12600, 12612), ("DDX11L1","pseudogene")) #Test start-in-bounds not in bounds of gene model entry

        #Weird values and new entries
        self.assertNotEqual(geneModel.identify_range('1', '+', -12700, -12800), ("DDX11L1","pseudogene")) #Handles negative loci
        self.assertEqual(geneModel.identify_range('MT','+', 15885, 15900), ("MT-CYB", "protein_coding")) #Test second chromosome entry
        self.assertNotEqual(geneModel.identify_range('MT','-',15885, 15900), ("MT-ND6", "protein_coding")) #Test alternate strand
        self.assertEqual(geneModel.identify_range('MT','+', 18001, 18005), ("MT-CYBz","protein_coding")) #Test additional entry

if __name__ == '__main__':
    unittest.main()

