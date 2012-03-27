'''
Created on Mar 27, 2012

@author: intern
'''

class Exon(object):
    def __init__(self, exon_id, query, target, length):
        self.id             = int(exon_id)
        self.query          = query
        self.target         = target
        self.length         = length
        
    def set_length (self, length):
        self.length         = length
        
    def set_intervals(self, q_start, q_end, t_start, t_end):
        self.q_start        = int(q_start)
        self.q_end          = int(q_end)
        self.t_start        = int(t_start)
        self.t_end          = int(t_end)
        
    def set_identity(self, match, alignmentLength):
        self.match          = int(match)
        self.alignmentLength= int(alignmentLength)
        
    def set_score(self, score):
        self.score          = int(score)    
        
    def __lt__ (self, other):
        if isinstance(other, Exon):
            significance_other = float(other.match) / other.length
            significance_self = float(self.match) / self.length
            if significance_other > significance_self:
                return True
            else:
                return False
        