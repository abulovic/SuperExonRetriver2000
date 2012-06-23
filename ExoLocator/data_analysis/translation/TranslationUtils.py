'''
Created on Jun 21, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

def translate_ensembl_exons (ensembl_exons):
    
    cdna = ""
    #determine the frame
    first_exon = ensembl_exons[0]
    if first_exon.relative_start == 0:
        new_frame = first_exon.frame
    else:
        frame_complement = (first_exon.relative_start - first_exon.frame) % 3
        if frame_complement == 1:
            new_frame = 2
        elif frame_complement == 2:
            new_frame = 1
        else:
            new_frame = 0 
        
    #create the cdna
    for ee in ensembl_exons:
        new_sequence = ee.sequence[ee.relative_start:ee.relative_stop]
        cdna += new_sequence
    al_exon = AlignmentExonPiece("ensembl", -1, "", cdna)
    al_exon.set_frame(new_frame)
    return al_exon



def patch_sequences(alignment_exon, actual_exon_len):
    '''
    Pads the alignment sequences with leftovers from the previous
    exon, with the right amount of Ns in the beginning (for the initial gaps)
    and the right amount of Ns in the end of the sequence
    '''
    
    genomic_cdna = Seq("", IUPAC.ambiguous_dna)
    exon_cdna    = Seq("", IUPAC.ambiguous_dna)
    
    
    # decide what to append to the cDNA
    (translation_start, translation_end) = (alignment_exon.alignment_info["sbjct_start"], 
                                            alignment_exon.alignment_info["sbjct_end"])
    genomic_cdna +=  ( "N"* (translation_start - 1))
    exon_cdna    +=  ("N" * (translation_start - 1))
    
    genomic_cdna += alignment_exon.alignment_info["query_seq"]
    exon_cdna    += alignment_exon.alignment_info["sbjct_seq"]
    
    genomic_cdna += "N" * (actual_exon_len - translation_end)
    exon_cdna    += "N" * (actual_exon_len - translation_end)
    
    return (genomic_cdna, exon_cdna)


def process_insertion_in_exon(human_seq, spec_seq, ordinal):
    if len(spec_seq) in (1,2,4,5):
        exon_type = "frameshift"
    else:
        exon_type = "insertion"
        
    al_exon = AlignmentExonPiece(exon_type, ordinal, human_seq, spec_seq)
    return al_exon
        


def process_insertion_free_region(human_seq, spec_seq, frame, ordinal):
    
    number_of_ungapped_sequences = len(re.split("-+", str(spec_seq)))
    pattern_string = "([ATGCN]+)"
    for i in range (0, number_of_ungapped_sequences-1):
        pattern_string += "([-]+)([ATGCN]+)"
    pattern = re.compile (pattern_string)
    sequences = re.match (pattern, str(spec_seq))
    
    al_exons = []
    
    in_coding = True
    start, stop = 0,0
    coding_len = 0
    
    for seq in sequences.groups():
        stop += len(seq)
        
        species = seq
        human   = human_seq[start:stop]
        
        if in_coding:
            if coding_len == 0:
                new_frame = frame
            else:
                frame_status = abs(coding_len - frame) % 3
                if frame_status == 1:
                    new_frame = 2
                elif frame_status == 2:
                    new_frame = 1
                else:
                    new_frame = 0
                
            exon = AlignmentExonPiece("coding", ordinal, human, species)
            exon.set_frame(new_frame)
            al_exons.append(exon)
            
        in_coding = not in_coding
        start = stop
        coding_len += len(seq)
            
    return al_exons


def split_exon_seq (alignment_exon, coding_exon):
    
    genomic_dna = alignment_exon.alignment_info["query_seq"]
    exon_dna    = alignment_exon.alignment_info["sbjct_seq"]
    
    number_of_ungapped_sequences = len(re.split("-+", str(exon_dna)))
    pattern_string = "([ATGCN]+)"
    for i in range (0, number_of_ungapped_sequences-1):
        pattern_string += "([-]+)([ATGCN]+)"
    pattern = re.compile (pattern_string)
    sequences = re.match (pattern, str(exon_dna))
    
    all_alignment_exons = []
    
    coding_len = 0
    start, stop = 0,0
    ordinal = 0
    in_coding = True
    
    # determine starting frame
    if alignment_exon.alignment_info["sbjct_start"] == 1:
        frame = coding_exon.frame
    elif alignment_exon.alignment_info["sbjct_start"] == 2:
        if coding_exon.frame == 0:
            frame = 2
        elif coding_exon.frame == 1:
            frame = 0
        else:
            frame = 1
    else:
        frame_status = abs(alignment_exon.alignment_info["sbjct_start"] - 1 - coding_exon.frame) % 3
        if frame_status == 2:
            frame = 1
        elif frame_status == 1:
            frame = 2
        else:
            frame = 0
 
    for seq in sequences.groups ():
        
        stop += len(seq)
        
        # determine the working sequences
        human_seq = seq
        spec_seq  = genomic_dna[start:stop]
        
        if not in_coding:
            exon = process_insertion_in_exon (human_seq, spec_seq, ordinal)
            all_alignment_exons.append(exon)

        # determine frame
        if in_coding:
            if coding_len == 0:
                new_frame = frame
            else:
                frame_status = abs(coding_len - frame) % 3
                if frame_status == 1:
                    new_frame = 2
                elif frame_status == 2:
                    new_frame = 1
                else:
                    new_frame = 0
            
        if in_coding:
            al_exons = process_insertion_free_region (human_seq, spec_seq, new_frame, ordinal)
            all_alignment_exons.extend(al_exons)
        
        
        # update the status' and lengths
        start = stop
        if in_coding:
            coding_len += len(seq)
        in_coding = not in_coding
        ordinal += 1
        
    return all_alignment_exons
        
    
    
    
    



def translate_alignment_exons(ref_protein_id, reference_species, alignment_exon):
    
    # load the reference exons
    exon_container  = ExonContainer.Instance()
    reference_exons = exon_container.get ((ref_protein_id, reference_species, "ensembl"))
    
    # get just the coding exons
    for coding_exon in reference_exons.get_coding_exons():
        if coding_exon.exon_id == alignment_exon.ref_exon_id:
            break
    coding_seq_len = len(coding_exon.sequence)
    
    
    (genomic_cdna, exon_cdna) = patch_sequences (alignment_exon, coding_seq_len)
    
    return (genomic_cdna, exon_cdna)
    
    
class BestExonAlignment (object):
    def __init__ (self, ref_exon_id, sw_gene_exon = None, ensembl_exons = None, status = None):
        self.ref_exon_id = ref_exon_id
        self.sw_gene_exon = sw_gene_exon
        self.ensembl_exons = ensembl_exons
        self.status = status
        
    def set_sw_gene_exons (self, sw_gene_exon):
        self.sw_gene_exon = sw_gene_exon
        
    def set_ensembl_exons (self, ensembl_exons):
        self.ensembl_exons= ensembl_exons
        
    def set_status (self, status):
        if status not in ("ensembl", "sw_gene", "both"):
            raise ValueError ("Status must be ensembl, sw_gene or both")
        self.status = status
        
    def set_scores (self, sw_gene_score, ensembl_score):
        self.sw_gene_score= sw_gene_score
        self.ensembl_score = ensembl_score
        
    def set_ref_sequences (self, gene_ref_seq, cdna_ref_seq):
        self.gene_ref_seq = gene_ref_seq
        self.cdna_ref_seq = cdna_ref_seq
        

class AlignmentExonPiece (object):
    '''
    Class representing one alignment piece. It is a part 
    of one exon alignment. It can be of various types:
     - coding: means that it corresponds to the coding part of the reference species DNA
     - insertion: meaning that it is an insertion in the coding ref species DNA
     - frameshift: meaning that is an insertion, but of length 1,2,4 or 5 base pairs
    '''
    def __init__ (self, exon_type, ordinal, ref_seq, spec_seq):
        '''
        @param exon_type: coding / insertion / frameshift
        @param ref_seq: nucleotide sequence in the reference species
        @param spec_seq: nucleotide sequence in the species
        @param ordinal: ordinal in the one exon alignment
        '''
        self.type       = exon_type
        self.ordinal    = ordinal
        self.ref_seq    = ref_seq
        self.spec_seq   = spec_seq
        
    def set_location (self, start, stop):
        '''
        Set location relative to the alignment
        '''
        self.start = start
        self.stop = stop
        
    def set_frame(self, frame):
        '''
        Set the translation frame (0,1,2)
        '''
        self.frame = frame
        
    
        








