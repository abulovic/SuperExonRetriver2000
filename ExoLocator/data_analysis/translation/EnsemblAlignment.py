'''
Created on Jun 24, 2012

@author: anana
'''
from data_analysis.containers.ProteinContainer import ProteinContainer
from Bio.Seq import Seq
import re
from data_analysis.translation.TranslationUtils import translate_ensembl_exons
from Bio.Alphabet import IUPAC


class EnsemblAlignment (object):
    
    def __init__ (self, ref_protein_id, ref_exon, alignment_exon, ensembl_exons):
        
        self.ref_protein_id = ref_protein_id
        self.ref_exon = ref_exon
        self.alignment_exon = alignment_exon
        self.ensembl_exons = ensembl_exons
        self.set_protein_sequences ()
        
    def set_protein_sequences (self):
        

        pc = ProteinContainer.Instance()
        
        ref_protein = pc.get(self.ref_protein_id)
        ref_protein_seq = ref_protein.get_sequence_record().seq
        partial_ref_seq = Seq(self.alignment_exon.alignment_info["sbjct_seq"].replace("-",""), IUPAC.ambiguous_dna)
        
        
        complete_protein_exon_seq = self.ref_exon.sequence [self.ref_exon.frame:].translate()
        if str(complete_protein_exon_seq).endswith("*"):
            complete_protein_exon_seq = complete_protein_exon_seq[0:len(complete_protein_exon_seq)-1]

        exon_prot_start = str(ref_protein_seq).find(str(complete_protein_exon_seq))
        exon_prot_stop = exon_prot_start + len(complete_protein_exon_seq)
        
        for frame in (0,1,2):
            partial_protein_ref_seq = partial_ref_seq[frame:].translate()
            if str(partial_protein_ref_seq).endswith("*"):
                partial_protein_ref_seq = partial_protein_ref_seq[0:len(partial_protein_ref_seq)-1]
            found = False
            if str(complete_protein_exon_seq).find (str(partial_protein_ref_seq)) != -1:
                for a in list(re.finditer(str(partial_protein_ref_seq), str(ref_protein_seq))): 
                    if a.start() >= exon_prot_start and a.end() <= exon_prot_stop:
                        self.ref_protein_seq = partial_protein_ref_seq
                        self.ref_protein_start = a.start()
                        self.ref_protein_stop = a.end()
                        found = True
                        break
            if found:
                break
                    
        self.spec_protein_seq = translate_ensembl_exons(self.ensembl_exons)