# module name: BedEntry
# main program: samplespecificdbgenerator

__author__ = "anthony"
__date__ = "$Oct 27, 2015 5:22:01 PM$"

from Bio.Seq import reverse_complement, translate
import sys

class BedEntry( object ):
    def __init__(self, line):
        self.line = line
        try:
            (chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts,seq) = line.split('\t')[0:13]
            self.chrom = chrom
            self.chromStart = int(chromStart)
            self.chromEnd = int(chromEnd)
            self.name = name
            self.score = int(score)
            self.strand = strand
            self.thickStart = int(thickStart)
            self.thickEnd = int(thickEnd)
            self.itemRgb = itemRgb
            self.blockCount = int(blockCount)
            self.blockSizes = [int(x) for x in blockSizes.split(',')]
            self.blockStarts = [int(x) for x in blockStarts.split(',')]
            self.seq = seq
        except Exception, e:
            print >> sys.stderr, "Unable to read Bed entry: %s" % e
            exit(1)

#    def get_splice_junctions(self):
#        splice_juncs = []
#        for i in range(self.blockCount  - 1):
#            splice_junc = "%s:%d_%d" % (self.chrom, self.chromStart + self.blockSizes[i], self.chromStart + self.blockStarts[i+1])
#            splice_juncs.append(splice_junc)
#        return splice_juncs

    def get_exon_seqs(self):
        exons = []
        for i in range(self.blockCount):
            # splice_junc = "%s:%d_%d" % (self.chrom, self.chromStart + self.blockSizes[i], self.chromStart + self.blockStarts[i+1])
            exons.append(self.seq[self.blockStarts[i] : self.blockStarts[i] + self.blockSizes[i]])
        if self.strand == '-':  #reverse complement
            exons.reverse()
            for i, s in enumerate(exons):
                exons[i] = reverse_complement(s)
        return exons
    
    def get_spliced_seq(self):
        return ''.join(self.get_exon_seqs()).strip()

    def get_translation(self, sequence):
        translation = None
        seq = sequence.strip()
        if seq and len(seq) >= 3: translation = translate(seq)
        return translation

    ## [[start,end,seq],[start,end,seq],[start,end,seq]]
    ## filter: ignore translation if stop codon in first exon after ignore_left_bp
    def get_filtered_translations(self):
        translations = [None, None, None]
        seq = self.get_spliced_seq()
        block_sum = sum(self.blockSizes)
        exon_sizes = self.blockSizes
        if self.strand == '-': exon_sizes.reverse()
        splice_sites = [sum(exon_sizes[:x]) / 3 for x in range(1, len(exon_sizes))] #exon_sizes can be assumed to be length 2
        junc = splice_sites[0] if len(splice_sites) > 0 else exon_sizes[0]
        if seq:
            for i in range(3):
                stop_translation = len(seq) - ((len(seq) - i) % 3)
                translation = self.get_translation(sequence = seq[i:stop_translation])
                if translation:
                    hasEx1StopCodon = translation.rfind('*', 0, junc)
                    if hasEx1StopCodon >= 0: continue
                    tstart = 0
                    stop = translation.find('*',junc)
                    tstop = stop if stop >= 0 else len(translation)
                    trimmed = translation[tstart:tstop]
                    #get genomic locations for start and end
                    offset = (block_sum - i) % 3
                    if self.strand == '+':
                        chromStart = self.chromStart + i + (tstart * 3)
                        chromEnd = self.chromEnd - offset - (len(translation) - tstop) * 3
                    else:
                        chromStart = self.chromStart + offset + (len(translation) - tstop) * 3
                        chromEnd = self.chromEnd - i - (tstart * 3)
                    translations[i] = [chromStart, chromEnd, trimmed]
        return translations

    def get_seq_id(self, seqtype='unk:unk', reference='', frame=None):
        ## Ensembl fasta ID format
        # >ID SEQTYPE:STATUS LOCATION GENE TRANSCRIPT
        # >ENSP00000328693 pep:splice chromosome:NCBI35:1:904515:910768:1 gene:ENSG00000158815:transcript:ENST00000328693 gene_biotype:protein_coding transcript_biotype:protein_coding
        frame_name, chromStart, chromEnd = '', self.chromStart, self.chromEnd
        strand = 1 if self.strand == '+' else -1
        if frame != None:
            block_sum = sum(self.blockSizes)
            offset = (block_sum - frame) % 3
            frame_name = '_' + str(frame + 1)
            if self.strand == '+':
                chromStart += frame
                chromEnd -= offset
            else:
                chromStart += offset
                chromEnd -= frame
        location = "chromosome:%s:%s:%s:%s:%s" % (reference, self.chrom, chromStart, chromEnd, strand)
        seq_id = "%s%s %s %s" % (self.name, frame_name, seqtype, location)
        return seq_id

    def get_line(self, start_offset = 0, end_offset = 0):
        if start_offset or end_offset:
            s_offset = start_offset if start_offset else 0
            e_offset = end_offset if end_offset else 0
            if s_offset > self.chromStart:
                s_offset = self.chromStart
            chrStart = self.chromStart - s_offset
            chrEnd = self.chromEnd + e_offset
            blkSizes = self.blockSizes
            blkSizes[0] += s_offset
            blkSizes[-1] += e_offset
            blkStarts = self.blockStarts
            for i in range(1,self.blockCount):
                blkStarts[i] += s_offset
            items = [str(x) for x in [self.chrom,chrStart,chrEnd,self.name,self.score,self.strand,self.thickStart,self.thickEnd,self.itemRgb,self.blockCount,','.join([str(x) for x in blkSizes]),','.join([str(x) for x in blkStarts])]]
            return '\t'.join(items) + '\n'
        return self.line