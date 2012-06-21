'''
Created on Jun 19, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.utilities.generate_structure import fill_all_containers

class BestProteinProduct (object):
    
    def __init__ (self, ref_protein_id, species, reference_species):
        
        self.ref_protein_id = ref_protein_id
        self.species        = species
        self.ref_species    = reference_species
        
    def load_alignments (self):
        
        exon_container      = ExonContainer.Instance()
        self.ref_exons      = exon_container.get((self.ref_protein_id, self.ref_species, "ensembl"))
        self.gene_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_gene"))
        self.cDNA_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_exon"))
        self.ensembl_exons  = exon_container.get((self.ref_protein_id, self.species, "ensembl"))
        
    def decide_on_best_exons (self):
        
        coding_exons = self.ref_exons.get_coding_exons()
        self.best_exons = {}
        
        
        for ce in coding_exons:
            
            gene_score, cdna_score = 0,0
            if ce.exon_id in self.gene_exons.alignment_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
                if gene_al_exon.viability:
                    gene_score = gene_al_exon.alignment_info["score"]
                    
            if ce.exon_id in self.cDNA_exons.alignment_exons:
                cdna_al_exon = self.cDNA_exons.alignment_exons[ce.exon_id][0]
                if cdna_al_exon.viability:
                    cdna_score = cdna_al_exon.alignment_info["score"]
                    (start,stop) = (cdna_al_exon.alignment_info["query_start"],
                                    cdna_al_exon.alignment_info["query_end"])
                    
            # if there is no alignment        
            if not cdna_score and not gene_score:
                best_exon_alignment = None
            # if the alignments are of the same quality
            elif cdna_score and gene_score and cdna_score == gene_score:
                
                best_exon_alignment = BestExonAlignment(ce.exon_id, 
                                                        gene_al_exon,
                                                        self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop),
                                                        "both" )
                
                
            elif cdna_score > gene_score:
                
                best_exon_alignment = BestExonAlignment(ce.exon_id, 
                                                        gene_al_exon,
                                                        self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop),
                                                        "ensembl" )
                
            else:
                
                if not cdna_score:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, 
                                                        gene_al_exon,
                                                        None,
                                                        "sw_gene" )
                else:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, 
                                                        gene_al_exon,
                                                        self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop),
                                                        "sw_gene" )
             
            best_exon_alignment.set_scores(gene_score, cdna_score)   
            self.best_exons[ce.exon_id] = best_exon_alignment
            
            
    def translate_best_exons_to_protein (self):
        
        
        ref_coding_exons = self.ref_exons.get_coding_exons()
        for ref_coding_exon in ref_coding_exons:
            print ref_coding_exon.frame
            
        
    
class BestExonAlignment (object):
    def __init__ (self, ref_exon_id, sw_gene_exons = None, ensembl_exons = None, status = None):
        self.ref_exon_id = ref_exon_id
        self.sw_gene_exons = sw_gene_exons
        self.ensembl_exons = ensembl_exons
        self.status = status
        
    def set_sw_gene_exons (self, sw_gene_exons):
        self.sw_gene_exons = sw_gene_exons
        
    def set_ensembl_exons (self, ensembl_exons):
        self.ensembl_exons= ensembl_exons
        
    def set_status (self, status):
        if status not in ("ensembl", "sw_gene", "both"):
            raise ValueError ("Status must be ensembl, sw_gene or both")
        self.status = status
        
    def set_scores (self, sw_gene_score, ensembl_score):
        self.sw_gene_score= sw_gene_score
        self.ensembl_score = ensembl_score
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
    
if __name__ == '__main__':
    fill_all_containers(True)
    bpp = BestProteinProduct("ENSP00000382758", "Gorilla_gorilla", "Homo_sapiens")
    bpp.load_alignments()
    bpp.decide_on_best_exons()
    
    be = bpp.best_exons
    bpp.translate_best_exons_to_protein()
    print
    
        
    
