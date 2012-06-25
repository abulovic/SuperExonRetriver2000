'''
Created on Jun 24, 2012

@author: anana
'''

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
        
    def set_translations (self, ref_protein_seq, spec_protein_seq):
        '''
        Set the protein translation of the alignment piece.
        These alignments need not completely correspond to the DNA sequence, 
        because some of them might be extended on both sides to account 
        for the translation frame (so as not to lose AAs)
        '''
        self.ref_protein_seq = ref_protein_seq
        self.spec_protein_seq = spec_protein_seq
        
    def set_protein_locations (self, ref_protein_start, ref_protein_stop):
        '''
        Set the locations of the translated piece on the 
        referent species protein
        '''
        self.ref_protein_start = ref_protein_start
        self.ref_protein_stop = ref_protein_stop
        
    def set_alignment_locations (self, start, stop):
        '''
        Set locations of the alignment on the species gene DNA.
        '''
        self.alignment_start = start
        self.alignment_stop  = stop
        
    def set_genomic_locations (self, start, stop, sequence_id):
        '''
        Set the genomic location and the ID of the species DNA sequence.
        '''
        self.genomic_start = start
        self.genomic_stop  = stop
        self.sequence_id = sequence_id
        