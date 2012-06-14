'''
Created on Mar 27, 2012

@author: intern
'''

class Exon_translation:
    '''
    Auxiliary class for constructing the translation from exon to protein
    '''
    
    def __init__ (self, referent_exon, alignment_exon):
        self.id = alignment_exon.ordinal
        self.length = referent_exon.length
        self.query = alignment_exon.alignment_info["query_seq"]
        self.target = alignment_exon.alignment_info["sbjct_seq"]
        self.set_intervals(alignment_exon.alignment_info["query_start"], 
                           alignment_exon.alignment_info["query_end"], 
                           alignment_exon.alignment_info["sbjct_start"], 
                           alignment_exon.alignment_info["sbjct_end"])
        self.set_identity(alignment_exon.alignment_info["identities"], alignment_exon.alignment_info["length"])
        self.set_viablity(alignment_exon.viability)
        
   
    def set_intervals(self, q_start, q_end, t_start, t_end):
        self.q_start    =   int(q_start)
        self.q_end      =   int(q_end)
        self.t_start    =   int(t_start)
        self.t_end      =   int(t_end)
    def set_identity(self, match, length):
        self.no_of_matches      =   int(match)
        self.alignment_length   =   int(length)
    def set_score(self, score):
        self.score              =   score
    def set_aligned_protein_ID (self, referece_species_ID):
        self.reference_species_ID = int(referece_species_ID)
    def set_viablity(self, viability):
        self.viability          =   viability
    def set_frame(self, frame):
        self.frame              =   frame
    def fitness(self):
        return float(self.alignment_length * self.no_of_matches) / self.length
        
    def __lt__ (self, other):
        if isinstance(other, Exon_translation):
            if other.fitness() > self.fitness():
                return True
            else:
                return False
        