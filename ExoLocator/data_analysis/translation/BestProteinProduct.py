'''
Created on Jun 19, 2012

@author: anana
'''

# BioPython imports
from Bio.Alphabet   import IUPAC
from Bio.Seq        import Seq

# data analysis imports
from data_analysis.containers.ExonContainer         import ExonContainer

from data_analysis.utilities.generate_structure     import fill_all_containers

from data_analysis.translation.EnsemblAlignment     import EnsemblAlignment
from data_analysis.translation.BestExonAlignment    import BestExonAlignment
from data_analysis.translation.SWGeneAlignment      import SWGeneAlignment


class BestProteinProduct (object):
    
    '''
    Class representing all the information necessary to reconstruct the
    best possible protein product. The best result can be offered
    either from Ensembl, from Exolocator or both (the same).
    '''
    
    def __init__ (self, ref_protein_id, species, reference_species):
        
        self.ref_protein_id = ref_protein_id
        self.species        = species
        self.ref_species    = reference_species
        self.load_alignments()
        
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
        self.gene_exons     = None
        try:
            self.gene_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_gene"))
        except KeyError:
            pass
        self.cDNA_exons     = None
        try:
            self.cDNA_exons     = exon_container.get((self.ref_protein_id, self.species, "sw_exon"))
        except KeyError:
            pass
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
            
            # refresh the auxiliary variables
            gene_score, cdna_score      = 0,0
            cdna_exons                  = None
            cdna_al_exon, gene_al_exon  = None, None
            gene_ref_seq, cdna_ref_seq  = "", ""
            ensembl_alignment           = None
            sw_gene_alignment           = None
            best_exon_alignment         = None
            
            ############## LOAD AVAILABLE DATA ############################################
            
            ##############    SW GENE DATA     ##########################################
            if self.gene_exons and ce.exon_id in self.gene_exons.alignment_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
                
                if gene_al_exon.viability:
                    gene_score = gene_al_exon.alignment_info["score"]
                    sw_gene_alignment = SWGeneAlignment(self.ref_protein_id, self.species, ce, gene_al_exon)

                else:
                    gene_al_exon = None
                        
            
            ##############    ENSEMBL DATA     ##########################################
            if self.ensembl_exons and self.cDNA_exons and ce.exon_id in self.cDNA_exons.alignment_exons:
                cdna_al_exon = self.cDNA_exons.alignment_exons[ce.exon_id][0]
                
                if cdna_al_exon.viability:
                    cdna_score = cdna_al_exon.alignment_info["score"]
                    (start,stop) = (cdna_al_exon.alignment_info["query_start"],
                                    cdna_al_exon.alignment_info["query_end"])
                    cdna_exons = self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop)
                    ensembl_alignment = EnsemblAlignment(self.ref_protein_id, ce, cdna_al_exon, cdna_exons)
                    
                else:
                    cdna_al_exon = None
                    
            
            ############## DETERMINE STATUS (ENSEMBL / SW_GENE / BOTH) ####################
                    
            # if there is no alignment        
            if not cdna_score and not gene_score:
                best_exon_alignment = None
                
            # if the alignments are of the same quality
            elif cdna_score and gene_score and cdna_score == gene_score:         
                best_exon_alignment = BestExonAlignment(ce.exon_id, "both" , ensembl_alignment, sw_gene_alignment)
            
            # if the ensembl alignment is better than the sw gene alignment           
            elif cdna_score > gene_score:
                best_exon_alignment = BestExonAlignment(ce.exon_id, "ensembl", ensembl_alignment, sw_gene_alignment )
                
            else:
                if not cdna_score:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, "sw_gene", None, sw_gene_alignment )
                else:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, "sw_gene", ensembl_alignment, sw_gene_alignment )
            if best_exon_alignment: 
                best_exon_alignment.set_scores(gene_score, cdna_score)   
                
            self.best_exons[ce.exon_id] = best_exon_alignment


            
    def patch_interexon_AAS (self):
        
        '''
        Patch the protein translations of exons. 
        The first alignment pieces are patched if they are on the beginning of the exon
        and if the last alignment piece of the last exon is the end of that exon. 
        The same goes for the last alignment pieces.
        '''
        ######### THIS FUNCTION IS IN NEED OF REFACTORING -- BUT NO TIME :( ################
        
        if not self.gene_exons:
            return
        
        coding_exons = self.ref_exons.get_coding_exons()
        
        for i in range (0, len(coding_exons)):
            
            this_exon       = coding_exons[i]
            alignment_exon  = self.gene_exons.alignment_exons[this_exon.exon_id][0]
            
            best_alignment  = self.best_exons[this_exon.exon_id]
            if not best_alignment or not best_alignment.sw_gene_alignment:
                continue
            
            add_beg_ref, add_end_ref    = "", ""
            add_beg_spec, add_end_spec  = "", ""
            
            ##### ONLY PATCH THE FIRST ALIGNMENT PIECE IF IT'S ON THE BEGINNING OF   ###########
            ##### THE EXON AND THE PREVIOUS EXON LAST ALIGNMENT PIECE IS ON          ###########
            ##### THE END OF THE EXON                                                ###########
            
            if alignment_exon.alignment_info["sbjct_start"] == 1:
                if i-1 in range (0, len(coding_exons)):
                    previous_exon       = coding_exons[i-1]
                    previous_al_exon    = self.gene_exons.alignment_exons[previous_exon.exon_id][0]
                    
                    if previous_al_exon and previous_al_exon.alignment_info["sbjct_end"] == len(this_exon.sequence):
                        how_much_to_take = (3 - this_exon.frame) % 3
                        add_beg_ref      = previous_al_exon.alignment_info["sbjct_seq"][len(previous_al_exon.alignment_info["sbjct_seq"])-how_much_to_take:]
                        add_beg_spec     = previous_al_exon.alignment_info["query_seq"][len(previous_al_exon.alignment_info["query_seq"])-how_much_to_take:]
                        
            ##### ONLY PATCH THE LAST ALIGNMENT PIECE IF IT'S ON THE END OF THE EXON  ########### 
            ##### AND THE NEXT EXON FIRST ALIGNMENT PIECE IS ON THE BEGINNING         ###########
            
            if alignment_exon.alignment_info["sbjct_end"] == len(this_exon.sequence):
                if i+1 in range (0, len(coding_exons)):
                    next_exon       = coding_exons[i+1]
                    next_al_exon    = self.gene_exons.alignment_exons[next_exon.exon_id][0]
                    
                    if next_al_exon and next_al_exon.alignment_info["sbjct_start"] == 1:
                        last_al_piece    = best_alignment.sw_gene_alignment.alignment_pieces[-1]
                        how_much_to_take = (3 - (len(last_al_piece.ref_seq) - last_al_piece.frame)) % 3
                        add_end_ref      = next_al_exon.alignment_info["sbjct_seq"][0:how_much_to_take]
                        add_end_spec     = next_al_exon.alignment_info["query_seq"][0:how_much_to_take]
                        
            # DETERMINE THE REF PROTEIN START / END LOCATION UPDATES            
            if add_beg_ref:
                start_update = 1
            else:
                start_update = 0
            if add_end_ref:
                end_update = 1
            else:
                end_update = 0
               
            
            # DO THE ACTUAL PROTEIN SEQUENCE & LOCATION UPDATES
            if len(best_alignment.sw_gene_alignment.alignment_pieces) == 1:
                al_piece = best_alignment.sw_gene_alignment.alignment_pieces[0]
                
                new_ref_cdna = add_beg_ref + al_piece.ref_seq + add_end_ref
                new_spec_cdna = add_beg_spec + al_piece.spec_seq + add_end_spec
                
                if add_beg_ref:
                    translation_frame = 0
                else:
                    translation_frame = al_piece.frame

                new_ref_prot                        = Seq(new_ref_cdna[translation_frame:], IUPAC.ambiguous_dna).translate()
                al_piece.ref_protein_seq            = new_ref_prot
                al_piece.ref_protein_start         -= start_update
                al_piece.ref_protein_stop          += end_update
                
                
                new_spec_prot                       = Seq(new_spec_cdna[translation_frame:], IUPAC.ambiguous_dna).translate()
                al_piece.spec_protein_seq           = new_spec_prot
                
            else:
                
                first_al_piece                      = best_alignment.sw_gene_alignment.alignment_pieces[0]
                last_al_piece                       = best_alignment.sw_gene_alignment.alignment_pieces[-1]
                
                if add_beg_ref:
                    translation_frame = 0
                else:
                    translation_frame = first_al_piece.frame
                
                new_ref_cdna                        = add_beg_ref + first_al_piece.ref_seq
                new_spec_cdna                       = add_beg_spec + first_al_piece.spec_seq
                
                new_ref_prot                        = Seq(new_ref_cdna[translation_frame:], IUPAC.ambiguous_dna).translate()
                first_al_piece.ref_protein_seq      = new_ref_prot
                first_al_piece.ref_protein_start   -= start_update
                
                
                new_spec_prot                       = Seq(new_spec_cdna, IUPAC.ambiguous_dna).translate()
                first_al_piece.spec_protein_seq     = new_spec_prot
                
                new_ref_cdna                        = last_al_piece.ref_seq + add_end_ref
                new_ref_prot                        = Seq(new_ref_cdna[last_al_piece.frame:], IUPAC).translate()
                last_al_piece.ref_protein_seq       = new_ref_prot
                last_al_piece.ref_protein_end      += end_update
                
                new_spec_cdna                       = last_al_piece.spec_seq + add_end_spec
                new_spec_prot                       = Seq(new_spec_cdna[last_al_piece.frame:], IUPAC).translate()
                last_al_piece.spec_protein_seq      = new_spec_prot



if __name__ == '__main__':

    fill_all_containers(True)
    bpp = BestProteinProduct("ENSP00000382758", "Ailuropoda_melanoleuca", "Homo_sapiens")
    bpp.load_alignments()
    bpp.decide_on_best_exons()
    bpp.patch_interexon_AAS()
    
    be = bpp.best_exons
    bpp.translate_best_exons_to_protein()
    print
    
        
    
