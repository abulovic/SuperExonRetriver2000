'''
Created on Apr 25, 2012

@author: marioot
'''

class Exon(object):
    '''
    classdocs
    '''
    def __init__(self, exon_type):
        '''
        Constructor
        '''
        self.exon_type = exon_type
        self.reference_exons = {}
        
    def set_alignment_info_for_ref_exon (self, identities, positives, gaps, start, stop, sequence, ref_exon_len, reference_exon_id):
        alignment_info = {"identities"  : identities, 
                          "positives"   : positives, 
                          "gaps"        : gaps, 
                          "start"       : start, 
                          "stop"        : stop,
                          "ref_exon_len": ref_exon_len,
                          "sequence":sequence}
        if not reference_exon_id in self.reference_exons:
            self.reference_exons[reference_exon_id] = [alignment_info]
        else:
            self.reference_exons[reference_exon_id].append(alignment_info)
            
    def set_reference_exon (self, ref_exon):
        for alignment_info in self.reference_exons[ref_exon.exon_id]:
            alignment_info["ref_exon"] = ref_exon
        
        