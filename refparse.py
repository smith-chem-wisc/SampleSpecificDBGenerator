# module name: refparse
# main program: samplespecificdbgenerator

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:54:44 PM$"

import variantcalls

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None : HTML_NS, "xsi" : XSI_NS}
UP = '{'+HTML_NS+'}'

def aa_abbrev_dict():
## dictionary for Amino Acid Abbreviations
    aa_abbrev_dict = dict()
    aa_abbrev_dict['Phe'] = 'F'
    aa_abbrev_dict['Leu'] = 'L'
    aa_abbrev_dict['Ser'] = 'S'
    aa_abbrev_dict['Tyr'] = 'Y'
    aa_abbrev_dict['Cys'] = 'C'
    aa_abbrev_dict['Trp'] = 'W'
    aa_abbrev_dict['Pro'] = 'P'
    aa_abbrev_dict['His'] = 'H'
    aa_abbrev_dict['Gln'] = 'Q'
    aa_abbrev_dict['Arg'] = 'R'
    aa_abbrev_dict['Ile'] = 'I'
    aa_abbrev_dict['Met'] = 'M'
    aa_abbrev_dict['Thr'] = 'T'
    aa_abbrev_dict['Asn'] = 'N'
    aa_abbrev_dict['Lys'] = 'K'
    aa_abbrev_dict['Val'] = 'V'
    aa_abbrev_dict['Ala'] = 'A'
    aa_abbrev_dict['Asp'] = 'D'
    aa_abbrev_dict['Glu'] = 'E'
    aa_abbrev_dict['Gly'] = 'G'
    return aa_abbrev_dict

def condense_xml_entry(entry):
    for element in entry:
        if element.tag not in [UP+'protein',UP+'accession',UP+'name',UP+'gene',UP+'organism',UP+'proteinExistence',UP+'depth',UP+'sequence',UP+'feature',UP+'dbReference']:
            entry.remove(element)
        elif element.get('type') != 'modified residue' and element.tag == UP+'feature': entry.remove(element)
        elif element.get('type') != 'Ensembl' and element.tag == UP+'dbReference': entry.remove(element)
        elif element.tag == UP+'organism':
            for field in element:
                if field.tag != UP+'name': element.remove(field)
        elif element.tag == UP+'protein':
            for name in element:
                if name.tag != UP+'recommendedName': element.remove(name)
        else: continue        

#Get header and seq from protein fasta using transcript or protein ensembl accession
def get_protein_fasta_seq(id, protein_fasta): 
    i = 0
    for item in protein_fasta[0]:
        if item.find(id) >= 0: return protein_fasta[0][i], protein_fasta[1][i]
        i += 1
    return None, None
        
def read_fasta_to_xml(root, refFasta):
    seq = ""
    line = refFasta.readline().strip()
    while line != "":
        if line.startswith(">"):
            line = line.split()
            acc, seqtype, chromosome = line[0][1:], line[1], line[2]
            geneId, transcriptId, addedInfo = line[3].split(':')[1], line[4].split(':')[1], ' '.join(line[5:])
            line = refFasta.readline().strip()
            while not line.startswith(">"):
                if line == "": break
                seq += line
                line = refFasta.readline().strip()
            if seq.find('*') < 0: variantcalls.enter_seqvar(root, acc, seqtype, chromosome, addedInfo, '',  '', geneId, transcriptId, seq)
            seq = ""
          
#Reads the headers and sequences of a fasta file into RAM
def read_protein_fasta(protein_fasta):
    proteinFasta = ([],[]) #headers, #sequences
    line = protein_fasta.readline()
    while line.startswith ('#'):
        line = protein_fasta.readline()
    sequence = ""
    while line != "":
        if line.startswith(">"):
            proteinFasta[0].append(line)
            line = protein_fasta.readline()
            while not line.startswith(">"):
                if line == "": break
                sequence += line
                line = protein_fasta.readline()
            proteinFasta[1].append(sequence.replace('\n','').replace('\r',''))
            sequence = ""
    protein_fasta.close()
    return proteinFasta

#Returns information in an XML entry in fasta duplet: (header, sequence)
def xml_to_fasta(entry):
    header = ">"
    if entry.tag == UP+'copyright': return None
    if entry.get('dataset') == 'Ensembl':
        accession, name, geneInfo, chromosome, geneId, transcriptId, addedInfo = None, None, None, None, None, None, None
        accession = entry.find(UP+'accession').text
        name = entry.find(UP + 'name').text
        geneInfo = entry.find(UP+'gene')
        if geneInfo != None: 
            geneInfo.getiterator(UP+'name')
            for item in geneInfo:
                if item.get('type') == 'coords': chromosome = item.text
                if item.get('type') == 'primary':
                    geneId = 'gene:' + item.text
                    transcriptId = 'transcript:' + item.find(UP+'transcript').text
        addedInfo = entry.find(UP + 'protein').find(UP + 'recommendedName').find(UP + 'fullName').text
        score = entry.find(UP + 'proteinExistence').find(UP + 'depth')
        score = 'depth:' + score.get('reads') if score != None else ''
        headerInfo = [x for x in [accession, name, chromosome, geneId, transcriptId, addedInfo, score] if x != None] #remove None elements
        header += ' '.join(headerInfo)
    else:
        if entry.get('dataset') == 'Swiss-Prot': database = 'sp'
        else: database = 'tr'
        accession = entry.find(UP + 'accession').text
        name = entry.find(UP + 'name').text
        accession = '|'.join([database, accession, name])
        organism = entry.find(UP + 'organism')
        if organism != None: organism = 'OS=' + organism.find(UP+'name').text
        geneName = entry.find(UP+ 'gene')
        if geneName != None: geneName = 'GN=' + geneName.find(UP+'name').text
        headerInfo = [x for x in [accession, organism, geneName] if x != None] 
        header += ' '.join(headerInfo)
    return(header.replace('\n','').replace('\r',''),  entry.find(UP+'sequence').text.replace('\n','').replace('\r',''))
    
###GENE MODEL STRUCTURE###
class Annotation():
    def __init__(self, name, biotype):
        self.name = name
        self.biotype = biotype
            
class Chromosome():   
    def __init__(self, id, strand):
        self.id = id
        self.strand = strand
        self.annotation_types = {}
        self.ranges = []
        
    def update_last_range_index(self, i):
        self.lastRangeIndex = i
        
    def identify_range(self, start, end):
        for range in self.ranges:
            chromStart, chromEnd, annotation = range
            startInBounds = start <= chromEnd and start >= chromStart
            endInBounds = end <= chromEnd and end >= chromStart
            spanning = start < chromStart and end > chromEnd
            if startInBounds or endInBounds or spanning: return (annotation.name, annotation.biotype)
        return ('', 'novel')
        
    def enter_annotation(self, chromStart, chromEnd, name, biotype):
        if (name, biotype) not in self.annotation_types: self.annotation_types[(name, biotype)] = Annotation(name, biotype)
        self.ranges.append((chromStart, chromEnd, self.annotation_types[(name, biotype)]))
        
class GeneModel():
    def __init__(self):
        self.chromosomes = {}
        
    def get_chrom(self, id, strand):
        if (id, strand) not in self.chromosomes: self.chromosomes[(id, strand)] = Chromosome(id, strand)
        return self.chromosomes[(id, strand)]
    
    def identify_range(self, chrom, strand, start, stop):
        return self.chromosomes[(chrom, strand)].identify_range(start, stop)
    
    def new_entry(self, line):
        line = line.split('\t')
        chromStart, chromEnd, chrom, strand = int(line[3]), int(line[4]), line[0], line[6]
        attribute_list = line[8].split('; ')
        attributes = {} #This is what's consistent across Ensembl gene models...
        for item in attribute_list:
            item = item.split(' ')
            attributes[item[0]] = item[1][1:-1]
        if 'transcript_biotype' not in attributes: attributes['transcript_biotype'] = line[1] #canonical gtf
        
        chromosome = self.get_chrom(chrom, strand)
        chromosome.enter_annotation(chromStart, chromEnd, attributes['gene_name'], attributes['gene_biotype'])