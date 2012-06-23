'''
Created on Jun 19, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from data_analysis.utilities.generate_structure import fill_all_containers
from data_analysis.translation.TranslationUtils import translate_ensembl_exons,\
    translate_alignment_exons, split_exon_seq, BestExonAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from utilities.DirectoryCrawler import DirectoryCrawler


def print_al_exons(al_exons, output_file):
    output_file.write("EXOLOCATOR:\n")
    for al_exon in al_exons:
        if al_exon.type == "coding":
            s1 = Seq(al_exon.ref_seq, IUPAC.ambiguous_dna)
            s2 = Seq(al_exon.spec_seq, IUPAC.ambiguous_dna)
            print "\tHUMAN:   ", s1[al_exon.frame:].translate()
            
            output_file.write ("\tHUMAN:       %s\n\tSPECIES:     %s\n" % (s1[al_exon.frame:].translate(), s2[al_exon.frame:].translate()))
            
            print "\tSPECIES: ", s2[al_exon.frame:].translate()
            


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
            gene_ref_seq, cdna_ref_seq = "", ""
            
            if not self.ensembl_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
            
            if ce.exon_id in self.gene_exons.alignment_exons:
                gene_al_exon = self.gene_exons.alignment_exons[ce.exon_id][0]
                if gene_al_exon.viability:
                    gene_score = gene_al_exon.alignment_info["score"]
                    gene_ref_seq = gene_al_exon.alignment_info["sbjct_seq"]
            
            # check cdna only if there exist ensembl exons for this species        
            if ce.exon_id in self.cDNA_exons.alignment_exons and self.ensembl_exons:
                cdna_al_exon = self.cDNA_exons.alignment_exons[ce.exon_id][0]
                if cdna_al_exon.viability:
                    cdna_score = cdna_al_exon.alignment_info["score"]
                    (start,stop) = (cdna_al_exon.alignment_info["query_start"],
                                    cdna_al_exon.alignment_info["query_end"])
                    cdna_exons = self.ensembl_exons.get_exon_ids_from_ccDNA_locations(start, stop)
                    cdna_ref_seq = cdna_al_exon.alignment_info["sbjct_seq"]
                    
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
                output_file.write ("\n\tTRANSLATION: %s\n" % str(al_exon.spec_seq[al_exon.frame:].translate()))
            if bea and bea.sw_gene_exon:
                (cdna_genomic, cdna_exon) = translate_alignment_exons(self.ref_protein_id, 
                                                                     self.ref_species, 
                                                                     bea.sw_gene_exon)
                al_exons = split_exon_seq(bea.sw_gene_exon, ref_coding_exon )
                print_al_exons(al_exons, output_file)
            print
            output_file.write ("\n")
            
        output_file.close()
            #species += cdna_genomic
            #human   += cdna_exon
   

if __name__ == '__main__':
    fill_all_containers(True)
    bpp = BestProteinProduct("ENSP00000275072", "Spermophilus_tridecemlineatus", "Homo_sapiens")
    bpp.load_alignments()
    bpp.decide_on_best_exons()
    
    be = bpp.best_exons
    bpp.translate_best_exons_to_protein()
    print
    
        
    
