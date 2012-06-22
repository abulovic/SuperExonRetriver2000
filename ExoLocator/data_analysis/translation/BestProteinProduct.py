'''
Created on Jun 19, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.utilities.generate_structure import fill_all_containers
from data_analysis.translation.TranslationUtils import translate_ensembl_exons,\
    translate_alignment_exons, split_exon_seq
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def print_al_exons(al_exons):
    for al_exon in al_exons:
        if al_exon.type == "coding":
            s1 = Seq(al_exon.ref_seq, IUPAC.ambiguous_dna)
            s2 = Seq(al_exon.spec_seq, IUPAC.ambiguous_dna)
            print "HUMAN:   ", s1[al_exon.frame:].translate()
            
            print "SPECIES: ", s2[al_exon.frame:].translate()


class BestProteinProduct (object):
    
    def __init__ (self, ref_protein_id, species, reference_species):
        
        self.ref_protein_id = ref_protein_id
        self.species        = species
        self.ref_species    = reference_species
        
    def load_alignments (self):
        '''
        Load all the necessary data, meaning:
         - reference exons
         - SW (gene to exon) alignment exons
         - SW (cdna to exon) alignment exons
         - ensembl species exons, if such exist
        '''
        
        exon_container      = ExonContainer.Instance()
        self.ref_exons      = exon_container.get((self.ref_protein_id, self.ref_species, "ensembl"))
        self.gene_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_gene"))
        self.cDNA_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_exon"))
        self.ensembl_exons  = None
        try:
            self.ensembl_exons  = exon_container.get((self.ref_protein_id, self.species, "ensembl"))
        except KeyError:
            pass
        
    def decide_on_best_exons (self):
        '''
        Map the best possible exon from the species to the 
        reference exon id. Best exons can come from the
        SW exon to gene alignment or from the Ensembl annotated exons
        '''
        
        coding_exons = self.ref_exons.get_coding_exons()
        self.best_exons = {}
         
        for ce in coding_exons:
            
            gene_score, cdna_score = 0,0
            cdna_exons = None
            cdna_al_exon, gene_al_exon = None, None
            
            if not self.ensembl_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
            
            if ce.exon_id in self.gene_exons.alignment_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
                if gene_al_exon.viability:
                    gene_score = gene_al_exon.alignment_info["score"]
            
            # check cdna only if there exist ensembl exons for this species        
            if ce.exon_id in self.cDNA_exons.alignment_exons and self.ensembl_exons:
                cdna_al_exon = self.cDNA_exons.alignment_exons[ce.exon_id][0]
                if cdna_al_exon.viability:
                    cdna_score = cdna_al_exon.alignment_info["score"]
                    (start,stop) = (cdna_al_exon.alignment_info["query_start"],
                                    cdna_al_exon.alignment_info["query_end"])
                    cdna_exons = self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop)
                    
            # if there is no alignment        
            if not cdna_score and not gene_score:
                best_exon_alignment = None
            # if the alignments are of the same quality
            elif cdna_score and gene_score and cdna_score == gene_score:             
                best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "both" )
                       
            elif cdna_score > gene_score:
                best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "ensembl" )
                
            else:
                if not cdna_score:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, None, "sw_gene" )
                else:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "sw_gene" )
            if best_exon_alignment: 
                best_exon_alignment.set_scores(gene_score, cdna_score)   
            self.best_exons[ce.exon_id] = best_exon_alignment
            
            
    def translate_best_exons_to_protein (self):
        
        species = Seq("", IUPAC.ambiguous_dna)
        human   = Seq("", IUPAC.ambiguous_dna)
        
        ref_coding_exons = self.ref_exons.get_coding_exons()
        for ref_coding_exon in ref_coding_exons:
            bea = self.best_exons[ref_coding_exon.exon_id]
            (extra_genomic, extra_exon) = ("", "")
            print "REF:     ", ref_coding_exon.sequence[ref_coding_exon.frame:].translate()
            if bea:
                print bea.status
            if bea and bea.ensembl_exons:
                protein = translate_ensembl_exons(bea.ensembl_exons)
                print "ENSEMBL: ", protein
            if bea and bea.sw_gene_exon:
                (cdna_genomic, cdna_exon) = translate_alignment_exons(self.ref_protein_id, 
                                                                     self.ref_species, 
                                                                     bea.sw_gene_exon)
                al_exons = split_exon_seq(bea.sw_gene_exon, ref_coding_exon )
                print_al_exons(al_exons)
            print
            #species += cdna_genomic
            #human   += cdna_exon
            
        
  
    
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
        
        

if __name__ == '__main__':
    fill_all_containers(True)
    bpp = BestProteinProduct("ENSP00000275072", "Otolemur_garnettii", "Homo_sapiens")
    bpp.load_alignments()
    bpp.decide_on_best_exons()
    
    be = bpp.best_exons
    bpp.translate_best_exons_to_protein()
    print
    
        
    
