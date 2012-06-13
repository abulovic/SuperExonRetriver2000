'''
Created on Apr 25, 2012

@author: marioot
'''

class Exon(object):
    '''
    Class representing one alignment exon. Take notice - alignment exon 
    does not need to correspont fully to the ensembl exon, it might be
    just a part of a whole exon.
    '''
    def __init__(self, exon_type, ref_exon_id, ref_protein_id, species):
        '''
        @param ref_exon_id: Ensembl ID of the exon which was aligned to the
        species genomic DNA resulting in this exon
        '''
        self.exon_type = exon_type
        self.ref_exon_id = ref_exon_id
        self.ref_protein_id = ref_protein_id
        self.species        = species
        
    def set_alignment_info (self, identities, positives, gaps, sbjct_start, sbjct_end, query_start, query_end, length, sbjct_seq, query_seq):
        self.alignment_info = {"identities"  : int(identities), 
                               "positives"   : int(positives), 
                               "gaps"        : int(gaps), 
                               "sbjct_start" : int(sbjct_start), 
                               "sbjct_end"   : int(sbjct_end),
                               "query_start" : int(query_start),
                               "query_end"   : int(query_end),
                               "length"      : int(length),
                               "sbjct_seq"   : sbjct_seq,
                               "query_seq"   : query_seq}
        
    def get_fitness (self):
        '''
        Measure of quality of an exon
        Quality is proportional to alignment length and number of identities, scaled by exon length
        '''
        return float(self.alignment_info["length"] * self.alignment_info["identities"]) / (self.alignment_info["query_end"]-self.alignment_info["query_start"])
    
    def set_viability (self, viability):
        '''
        Exon is viable if it is in the correct order with respect to surrounding exons
        '''
        self.viability = viability
        
    def set_ordinal (self, ordinal, alignment_ordinal):
        '''
        Set exon ordinal
        '''
        self.ordinal = ordinal
        self.alignment_ordinal = alignment_ordinal
        