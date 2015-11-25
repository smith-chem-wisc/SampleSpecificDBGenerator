# test module: test_novelsplices.py
# main program: samplespecificdbgenerator.py

import unittest
import novelsplices
import refparse
import BedEntry
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

class  Test_enter_seqvar(unittest.TestCase):
    def setUp(self):
        self.root = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
        self.db = et.ElementTree(self.root)

    def test_enter_seqvar(self):
        novelsplices.enter_seqvar(self.root, ">accession", "peptide:seqtype", "chromosome", "full name", "loci", "score", "geneId", "transcriptId", "seq\nuence")
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
        
class  Test_generate_tryptic_peptides(unittest.TestCase):
    def test_no_k_or_r(self):
        result = novelsplices.generate_tryptic_peps("AAA")
        self.assertEqual(result, ["AAA"])

    def test_k_nterminus(self):
        result = novelsplices.generate_tryptic_peps("KAAA")
        self.assertEqual(result, ["K", "AAA"])
        
    def test_k_cterminus(self):
        result = novelsplices.generate_tryptic_peps("AAAK")
        self.assertEqual(result, ["AAAK"])
        
    def test_r_nterminus(self):
        result = novelsplices.generate_tryptic_peps("RAAA")
        self.assertEqual(result, ["R", "AAA"])
        
    def test_r_cterminus(self):
        result = novelsplices.generate_tryptic_peps("AAAR")
        self.assertEqual(result, ["AAAR"])
        
    def test_internal_k_split(self):
        result = novelsplices.generate_tryptic_peps("KAAKAAAK")
        self.assertEqual(result, ["K", "AAK", "AAAK"])
        
    def test_internal_r_split(self):
        result = novelsplices.generate_tryptic_peps("RAARAAAR")
        self.assertEqual(result, ["R", "AAR", "AAAR"])
        
    def test_mix_split(self):
        result = novelsplices.generate_tryptic_peps("KARAAKAAARAAAAKAAAAA")
        self.assertEqual(result, ["K", "AR", "AAK", "AAAR", "AAAAK", "AAAAA"])
        
    def test_all_amino_acids(self):
        result = novelsplices.generate_tryptic_peps("FLSYCWPHQRIMTNKVADEG")
        self.assertEqual(result, ["FLSYCWPHQR", "IMTNK", "VADEG"])
        
    def test_non_standard_amino_acid(self):
        result = novelsplices.generate_tryptic_peps("X")
        self.assertEqual(result, ["X"])
        
class Test_update_tryptic_peptide_chromosome_indices(unittest.TestCase):
    def setUp(self):
        self.trypFragIndex = 10
        self.exon1Right = 20 
        self.exon2Left = 100
        self.peptide = "AAAAA" #5mer -> 15 nts; should end at 24
    
    #No boundary crossing
    def test_update_index_within_exon1(self):
        trypFragIndex = 10
        exon1Right = 20 
        exon2Left = 100
        peptide = "A" #[10,11,12],13
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 13)
        
    def test_update_index_2_before_exon1_right_boundary(self):
        trypFragIndex = 10
        exon1Right = 15 
        exon2Left = 100
        peptide = "A" #[10,11,12],13
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 13)
   
    def test_update_index_1_before_exon1_right_boundary(self):
        trypFragIndex = 10
        exon1Right = 14 
        exon2Left = 100
        peptide = "A" #[10,11,12],13
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 13)
        
    def test_update_index_up_to_exon1_right_boundary(self):
        trypFragIndex = 10
        exon1Right = 13 
        exon2Left = 100
        peptide = "A" #[10,11,12],13
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 13)
    
    #Crossing boundary
    def test_update_index_crossing_boundary(self):
        trypFragIndex = 10
        exon1Right = 13 
        exon2Left = 100
        peptide = "AA" #[10,11,12],[13,100,101],102
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 102)
        
    def test_update_index_1_more_than_boundary(self):
        trypFragIndex = 10
        exon1Right = 12 
        exon2Left = 100
        peptide = "A" #[10,11,12],100
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 100)
        
    def test_update_index_2_more_than_boundary(self):
        trypFragIndex = 10
        exon1Right = 11 
        exon2Left = 100
        peptide = "A" #[10,11,100],101
        trypFragIndex = novelsplices.update_tryp_index(trypFragIndex, exon1Right, exon2Left, peptide)
        self.assertEqual(trypFragIndex, 101)

if __name__ == '__main__':
    unittest.main()

