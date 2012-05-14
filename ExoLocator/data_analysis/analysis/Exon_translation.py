'''
Created on Mar 27, 2012

@author: intern
'''

class Exon_translation:
    def __init__(self, exon_id, length, query, target):
        self.id         =   int(exon_id)
        self.length     =   int(length)
        self.query      =   query
        self.target     =   target
        self.start_nuc  =   ""
        self.end_nuc    =   ""
        self.score      =   0
        self.no_of_matches      =   0
        self.alignment_length   =   0
        self.viability   =   True
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
        if isinstance(other, Exon):
            if other.fitness() > self.fitness():
                return True
            else:
                return False
        