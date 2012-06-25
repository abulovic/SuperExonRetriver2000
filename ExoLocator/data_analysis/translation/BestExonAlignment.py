'''
Created on Jun 24, 2012

@author: anana
'''
class BestExonAlignment (object):
    '''
    Class containing all the alignments for one exon, with a status
    describing which of the alignments is the best
    '''

    def __init__ (self, ref_exon_id, status = None, ensembl_alignment = None, sw_gene_alignment = None):
        self.ref_exon_id = ref_exon_id
        self.status         = status
        self.ensembl_alignment = ensembl_alignment
        self.sw_gene_alignment = sw_gene_alignment
        
    def set_sw_gene_exon (self, sw_gene_exon):
        self.sw_gene_exon   = sw_gene_exon
        
    def set_ensembl_exons (self, ensembl_exons):
        self.ensembl_exons  = ensembl_exons
        
    def set_status (self, status):
        '''
        @param status: ensembl / sw_gene / both
        '''
        if status not in ("ensembl", "sw_gene", "both"):
            raise ValueError ("Status must be ensembl, sw_gene or both")
        self.status = status
        
    def set_scores (self, sw_gene_score, ensembl_score):
        self.sw_gene_score = sw_gene_score
        self.ensembl_score = ensembl_score

        

