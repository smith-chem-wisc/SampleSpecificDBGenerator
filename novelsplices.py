# module name: novelsplices
# main program: samplespecificdatabases

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:55:01 PM$"

from lxml import etree as et
from BedEntry import BedEntry

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
    fullName.text = addedInfo

    gene = et.SubElement(entry, UP+'gene')
    geneName1 = et.SubElement(gene, UP+'name', type="coords")
    if loci: geneName1.text = chromosome + ":" + loci
    else: geneName1.text = chromosome

    organism = et.SubElement(entry, UP+'organism')
    organismName1 = et.SubElement(organism, UP+'name', type="scientific")
    organismName2 = et.SubElement(organism, UP+'name', type="common")
    organismName1.text = "Homo sapiens"
    organismName2.text = "Human"

    proteinExist = et.SubElement(entry, UP+'proteinExistence', type="evidence at transcript level")
    if score: et.SubElement(proteinExist, UP+'depth', reads=str(score))

    sequence = et.SubElement(entry, UP+'sequence', version="1", fragment="single")
    sequence.text = seq.replace('\n','').replace('\r','')

def generate_tryptic_peps(pep_seq):
    tryptic_peps = []
    k_fragments = pep_seq.split('K')
    for i, k_fragment in enumerate(k_fragments):
        if i != len(k_fragments) - 1: kr_fragments = (k_fragment + 'K').split('R')
        else: kr_fragments = k_fragment.split('R')
        for j, kr_fragment in enumerate(kr_fragments):
            if j != len(kr_fragments) - 1: tryptic_peps.append(kr_fragment + 'R')
            else: tryptic_peps.append(kr_fragment)
    while tryptic_peps and not tryptic_peps[-1]: tryptic_peps.pop() #remove empty elements at end of list
    return tryptic_peps

def update_tryp_index(trypIndex, exon1Right, exon2Left, peptide):
    trypIndex += len(peptide) * 3
    if trypIndex > exon1Right and trypIndex < exon2Left: #excedes the first exon 
        trypIndex = trypIndex - exon1Right + exon2Left - 1 #set to position in exon2
    return trypIndex

def translate_bed_line(root, geneModel, line, nsjDepthCut, minLength, refName, scoreName):
        entry = BedEntry(line)
        if entry.score <= nsjDepthCut: return # evaluate for depth cutoff
        translations = entry.get_filtered_translations()
        for i, tx in enumerate(translations):
            if tx:
                (chromStart, chromEnd, translation) = tx

                #Exon limit positions
                exon1Left = entry.chromStart
                exon2Left = entry.chromStart + entry.blockStarts[1]
                exon1Right = exon1Left + entry.blockSizes[0] - 1
                exon2Right = exon2Left + entry.blockSizes[1] - 1
                
                #Iterate through tryptic fragments and enter the one containing the junction into the database
                trypFragLeft, trypFragRight = chromStart, chromStart - 1
                trypticPeptides = generate_tryptic_peps(translation)
                for j, peptide in enumerate(trypticPeptides):
                    trypFragRight = update_tryp_index(trypFragRight, exon1Right, exon2Left, peptide)
                    
                    #Contains the junction, so look up type in gene model, e.g. exon1::novel;exon2:FOX1:protein_coding.
                    #Also, must not contain the arbitrary left end, which could be ragged and not tryptic. If it contains the right side, it must end with K or R.
                    if trypFragLeft < exon1Right and trypFragRight > exon1Right and len(peptide) >= minLength: 
                        exon1Type = geneModel.identify_range(entry.chrom, entry.strand, trypFragLeft, exon1Right)
                        exon2Type = geneModel.identify_range(entry.chrom, entry.strand, exon2Left, trypFragRight)

                        frame_name = '_%s' % (i + 1)
                        accession = entry.name + frame_name
                        chromosome = "chromosome:%s:%s:%s" % (refName, entry.chrom, entry.strand)
                        exons = "exons:%s:%s;%s:%s" % (exon1Type[0], exon1Type[1], exon2Type[0], exon2Type[1])
                        loci = "%s:%s;%s:%s" % (trypFragLeft, exon1Right, exon2Left, trypFragRight)
                        score = str(entry.score) if scoreName == 'depth' else ''
                        enter_seqvar(root, accession, 'pep:splice', chromosome, exons, loci, score, '', '', peptide)

                    trypFragLeft = update_tryp_index(trypFragLeft, exon1Right, exon2Left, peptide)
