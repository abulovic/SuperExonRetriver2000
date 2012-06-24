'''
Created on Jun 24, 2012

@author: intern
'''
class BestExonAlignment (object):
    '''
    Class containing all the alignments for one exon, with a status
    describing which of the alignments is the best
    '''
    def __init2__ (self, ref_exon_id, gene_alignment_exon, cdna_alignment_exon, ensembl_exons):
        pass
    def __init__ (self, ref_exon_id, sw_gene_exon = None, ensembl_exons = None, status = None, ensembl_alignment = None, sw_gene_alignment = None):
        self.ref_exon_id = ref_exon_id
        self.sw_gene_exon   = sw_gene_exon
        self.ensembl_exons  = ensembl_exons
        self.status         = status
        self.ensembl_alignment = ensembl_alignment
        self.sw_gene_alignment = sw_gene_alignment
        
    def set_sw_gene_exons (self, sw_gene_exon):
        self.sw_gene_exon   = sw_gene_exon
        
    def set_ensembl_exons (self, ensembl_exons):
        self.ensembl_exons  = ensembl_exons
        
    def set_status (self, status):
        if status not in ("ensembl", "sw_gene", "both"):
            raise ValueError ("Status must be ensembl, sw_gene or both")
        self.status = status
        
    def set_scores (self, sw_gene_score, ensembl_score):
        self.sw_gene_score = sw_gene_score
        self.ensembl_score = ensembl_score
        
    def set_ref_sequences (self, gene_ref_seq, cdna_ref_seq):
        self.gene_ref_seq = gene_ref_seq
        self.cdna_ref_seq = cdna_ref_seq
        
        

