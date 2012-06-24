'''
Created on Jun 19, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.utilities.generate_structure import fill_all_containers
from data_analysis.translation.TranslationUtils import translate_ensembl_exons, split_exon_seq
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from utilities.DirectoryCrawler import DirectoryCrawler
from data_analysis.translation.EnsemblAlignment import EnsemblAlignment
from data_analysis.translation.BestExonAlignment import BestExonAlignment
from data_analysis.translation.SWGeneAlignment import SWGeneAlignment


def print_al_exons(al_exons, output_file):
    output_file.write("EXOLOCATOR:\n")
    for al_exon in al_exons:
        if al_exon.type == "coding":
            s1 = Seq(al_exon.ref_seq, IUPAC.ambiguous_dna)
            s2 = Seq(al_exon.spec_seq, IUPAC.ambiguous_dna)

            print "\tHUMAN:   ", al_exon.ref_protein_seq
            print "\tprot: %d-%d, genome: %d-%d, %s" % (al_exon.ref_protein_start, al_exon.ref_protein_stop, al_exon.genomic_start, al_exon.genomic_stop, al_exon.sequence_id)
            #print "\tH(old):  ", Seq(al_exon.ref_seq[al_exon.frame:], IUPAC.ambiguous_dna).translate()
            output_file.write ("\tHUMAN:       %s\n\tSPECIES:     %s\n" % (s1[al_exon.frame:].translate(), s2[al_exon.frame:].translate()))
            
            print "\tSPECIES: ", al_exon.spec_protein_seq
            #print "\tS(old):  ", Seq(al_exon.spec_seq[al_exon.frame:], IUPAC.ambiguous_dna).translate()
            


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
            
            # if there is sw gene alignment, load all the data
            if self.gene_exons and ce.exon_id in self.gene_exons.alignment_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
                
                if gene_al_exon.viability:
                    gene_score = gene_al_exon.alignment_info["score"]
                    gene_ref_seq = gene_al_exon.alignment_info["sbjct_seq"]
                    sw_gene_alignment = SWGeneAlignment(self.ref_protein_id, self.species, ce, gene_al_exon)

                else:
                    gene_al_exon = None
                        
            
            # check cdna only if there exist ensembl exons for this species        
            if self.ensembl_exons and ce.exon_id in self.cDNA_exons.alignment_exons:
                cdna_al_exon = self.cDNA_exons.alignment_exons[ce.exon_id][0]
                
                if cdna_al_exon.viability:
                    cdna_score = cdna_al_exon.alignment_info["score"]
                    (start,stop) = (cdna_al_exon.alignment_info["query_start"],
                                    cdna_al_exon.alignment_info["query_end"])
                    cdna_exons = self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop)
                    cdna_ref_seq = cdna_al_exon.alignment_info["sbjct_seq"]
                    ensembl_alignment = EnsemblAlignment(self.ref_protein_id, ce, cdna_al_exon, cdna_exons)
                    
                else:
                    cdna_al_exon = None
                    
            
            ############## DETERMINE STATUS (ENSEMBL / SW_GENE / BOTH) ####################
                    
            # if there is no alignment        
            if not cdna_score and not gene_score:
                best_exon_alignment = None
            # if the alignments are of the same quality
            elif cdna_score and gene_score and cdna_score == gene_score:         
                best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "both" , ensembl_alignment, sw_gene_alignment)
                       
            elif cdna_score > gene_score:
                best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "ensembl", ensembl_alignment, sw_gene_alignment )
                
            else:
                if not cdna_score:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, None, "sw_gene", None, sw_gene_alignment )
                else:
                    best_exon_alignment = BestExonAlignment(ce.exon_id, gene_al_exon, cdna_exons, "sw_gene", ensembl_alignment, sw_gene_alignment )
            if best_exon_alignment: 
                best_exon_alignment.set_scores(gene_score, cdna_score)   
                best_exon_alignment.set_ref_sequences(gene_ref_seq, cdna_ref_seq)
                
            self.best_exons[ce.exon_id] = best_exon_alignment

            
            
    def translate_best_exons_to_protein (self):
        
        species = Seq("", IUPAC.ambiguous_dna)
        human   = Seq("", IUPAC.ambiguous_dna)
        
        dc = DirectoryCrawler()
        
        output_file_path = "%s/%s.fa" % (dc.get_mafft_path(self.ref_protein_id), self.species)
        output_file = open(output_file_path, "w")
        
        ref_coding_exons = self.ref_exons.get_coding_exons()
        for ref_coding_exon in ref_coding_exons:
            
            output_file.write ("REF_EXON: %s, FRAME: %d, ORDINAL: %d\n" % (ref_coding_exon.exon_id, ref_coding_exon.frame, ref_coding_exon.ordinal))
            output_file.write (str(ref_coding_exon.sequence[ref_coding_exon.frame:].translate()) + "\n")
            
            bea = self.best_exons[ref_coding_exon.exon_id]
            if bea:
                print "-------------", bea.ref_exon_id
                print bea.status
                output_file.write ("SOLUTION: %s\n" % bea.status)
            else:
                output_file.write ("SOLUTION: NO_SOLUTION\n")
                continue
            
            if bea and bea.ensembl_exons:
                al_exon = translate_ensembl_exons(bea.ensembl_exons)
                output_file.write("ENSEMBL:\n")
                output_file.write("\tEXONS_INCLUDED: ")
                for exon in bea.ensembl_exons:
                    output_file.write("%s:%d-%d, " % (exon.exon_id, exon.relative_start, exon.relative_stop))
                #output_file.write ("\n\tTRANSLATION: %s\n" % str(al_exon.spec_seq[al_exon.frame:].translate()))
                output_file.write ("\n\tTRANSLATION2:%s\n" % str(bea.ensembl_alignment.spec_protein_seq))
                #print "\n\tTRANSLATION: %s\n" % str(al_exon.spec_seq[al_exon.frame:].translate())
                print "\n\tHUMAN:       %s" % str(bea.ensembl_alignment.ref_protein_seq)
                print "\tTRANSLATION2:%s\n" % str(bea.ensembl_alignment.spec_protein_seq), bea.ensembl_alignment.ref_protein_start, bea.ensembl_alignment.ref_protein_stop
            if bea and bea.sw_gene_exon:
                al_exons = split_exon_seq(bea.sw_gene_exon, ref_coding_exon )
                if bea.sw_gene_exon:
                    print_al_exons(bea.sw_gene_alignment.alignment_pieces, output_file)
            print
            output_file.write ("\n")
            
        output_file.close()
            #species += cdna_genomic
            #human   += cdna_exon
            
    def patch_interexon_AAS (self):
        
        if not self.gene_exons:
            return
        
        coding_exons = self.ref_exons.get_coding_exons()
        for i in range (0, len(coding_exons)):
            
            this_exon = coding_exons[i]
            alignment_exon = self.gene_exons.alignment_exons[this_exon.exon_id][0]
            
            best_alignment = self.best_exons[this_exon.exon_id]
            if not best_alignment or not best_alignment.sw_gene_alignment:
                continue
            
            add_beg_ref, add_end_ref = "", ""
            add_beg_spec, add_end_spec = "", ""
            
            if alignment_exon.alignment_info["sbjct_start"] == 1:
                if i-1 in range (0, len(coding_exons)):
                    previous_exon = coding_exons[i-1]
                    previous_al_exon = self.gene_exons.alignment_exons[previous_exon.exon_id][0]
                    if previous_al_exon and previous_al_exon.alignment_info["sbjct_end"] == len(this_exon.sequence):
                        how_much_to_take = (3 - this_exon.frame) % 3
                        add_beg_ref  = previous_al_exon.alignment_info["sbjct_seq"][len(previous_al_exon.alignment_info["sbjct_seq"])-how_much_to_take:]
                        add_beg_spec = previous_al_exon.alignment_info["query_seq"][len(previous_al_exon.alignment_info["query_seq"])-how_much_to_take:]
            
            if alignment_exon.alignment_info["sbjct_end"] == len(this_exon.sequence):
                if i+1 in range (0, len(coding_exons)):
                    next_exon = coding_exons[i+1]
                    next_al_exon = self.gene_exons.alignment_exons[next_exon.exon_id][0]
                    if next_al_exon and next_al_exon.alignment_info["sbjct_start"] == 1:
                        last_al_piece = best_alignment.sw_gene_alignment.alignment_pieces[-1]
                        how_much_to_take = (3 - (len(last_al_piece.ref_seq) - last_al_piece.frame)) % 3
                        add_end_ref  = next_al_exon.alignment_info["sbjct_seq"][0:how_much_to_take]
                        add_end_spec = next_al_exon.alignment_info["query_seq"][0:how_much_to_take]
                        
            if add_beg_ref and add_end_ref:
                if len(best_alignment.sw_gene_alignment.alignment_pieces) == 1:
                    al_piece = best_alignment.sw_gene_alignment.alignment_pieces[0]
                    
                    new_ref_cdna = add_beg_ref + al_piece.ref_seq + add_end_ref
                    new_prot = Seq(new_ref_cdna, IUPAC.ambiguous_dna).translate()
                    al_piece.ref_protein_seq = new_prot
                    al_piece.ref_protein_start -= 1
                    al_piece.ref_protein_stop += 1
                    
                    new_spec_cdna = add_beg_spec + al_piece.spec_seq + add_end_spec
                    new_prot = Seq(new_spec_cdna, IUPAC.ambiguous_dna).translate()
                    al_piece.spec_protein_seq = new_prot
                    
                else:
                    first_al_piece = best_alignment.sw_gene_alignment.alignment_pieces[0]
                    last_al_piece = best_alignment.sw_gene_alignment.alignment_pieces[-1]
                    
                    new_ref_cdna = add_beg_ref + first_al_piece.ref_seq
                    new_prot = Seq(new_ref_cdna, IUPAC.ambiguous_dna).translate()
                    first_al_piece.ref_protein_seq = new_prot
                    first_al_piece.ref_protein_start -= 1
                    
                    new_spec_cdna = add_beg_spec + first_al_piece.spec_seq
                    new_prot = Seq(new_spec_cdna, IUPAC.ambiguous_dna).translate()
                    first_al_piece.spec_protein_seq = new_prot
                    
                    new_ref_cdna = last_al_piece.ref_seq + add_end_ref
                    new_prot = Seq(new_ref_cdna[last_al_piece.frame:], IUPAC).translate()
                    last_al_piece.ref_protein_seq = new_prot
                    last_al_piece.ref_protein_end += 1
                    
                    new_spec_cdna = last_al_piece.spec_seq + add_end_spec
                    new_prot = Seq(new_spec_cdna[last_al_piece.frame:], IUPAC).translate()
                    last_al_piece.spec_protein_seq = new_prot
                    
            elif add_beg_ref and not add_end_ref:
                al_piece = best_alignment.sw_gene_alignment.alignment_pieces[0]
                
                new_ref_cdna = add_beg_ref + al_piece.ref_seq
                new_prot = Seq(new_ref_cdna, IUPAC.ambiguous_dna).translate()
                al_piece.ref_protein_seq = new_prot
                al_piece.ref_protein_start -= 1
                
                new_spec_cdna = add_beg_spec + al_piece.spec_seq
                new_prot = Seq(new_spec_cdna, IUPAC.ambiguous_dna).translate()
                al_piece.spec_protein_seq = new_prot
                
            elif add_end_ref and not add_beg_ref:
                al_piece = best_alignment.sw_gene_alignment.alignment_pieces[-1]
                
                new_ref_cdna = al_piece.ref_seq + add_end_ref
                new_prot = Seq(new_ref_cdna[al_piece.frame:], IUPAC.ambiguous_dna).translate()
                al_piece.ref_protein_seq = new_prot
                al_piece.ref_protein_stop += 1
                
                new_spec_cdna = al_piece.spec_seq + add_end_spec
                new_prot = Seq(new_spec_cdna[al_piece.frame:], IUPAC.ambiguous_dna).translate()
                al_piece.spec_protein_seq = new_prot
                    
                    


if __name__ == '__main__':

    fill_all_containers(True)
    bpp = BestProteinProduct("ENSP00000382758", "Ailuropoda_melanoleuca", "Homo_sapiens")
    bpp.load_alignments()
    bpp.decide_on_best_exons()
    bpp.patch_interexon_AAS()
    
    be = bpp.best_exons
    bpp.translate_best_exons_to_protein()
    print
    
        
    
