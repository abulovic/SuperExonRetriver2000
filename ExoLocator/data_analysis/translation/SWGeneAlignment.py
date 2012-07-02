'''
Created on Jun 24, 2012

@author: anana
'''
# Python imports
import re

# utilities imports
from utilities.ConfigurationReader import ConfigurationReader

# data analysis imports
from data_analysis.translation.TranslationUtils import process_exon_alignment, set_protein_sequences

from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer


class SWGeneAlignment (object):
    
    def __init__ (self, ref_protein_id, species, ref_exon, alignment_exon):   
        
        self.ref_protein_id = ref_protein_id
        self.species        = species
        self.ref_exon       = ref_exon
        self.alignment_exon = alignment_exon
        self.alignment_pieces    = process_exon_alignment(self.alignment_exon, self.ref_exon)
        self.set_protein_sequences()
        self.determine_absolute_coordinates ()
    
            
    def set_protein_sequences (self):
        '''
        Loads the alignment pieces and sets their
        translations to protein.
        '''
        
        pc = ProteinContainer.Instance()
        
        ref_protein             = pc.get(self.ref_protein_id)
        ref_protein_seq         = ref_protein.get_sequence_record().seq
        ref_exon_translation    = self.ref_exon.sequence[self.ref_exon.frame:].translate()
        
        # remove the stop codon from the last position
        if str(ref_exon_translation).endswith("*"):
            ref_exon_translation = ref_exon_translation[0:len(ref_exon_translation)-1]
        
        self.alignment_pieces    = process_exon_alignment(self.alignment_exon, self.ref_exon)
        
        self.alignment_pieces    = set_protein_sequences (self.alignment_pieces)
        
        # find the locations of the exon translation in the protein
        exon_start = str(ref_protein_seq).find(str(ref_exon_translation))
        exon_stop = exon_start + len(ref_exon_translation)
   
        previous = None
   
        for al_piece in self.alignment_pieces:
            
            if al_piece.type == "coding":

                ref_protein_seq_piece = str(al_piece.ref_protein_seq)
                if ref_protein_seq_piece.endswith("*"):
                    ref_protein_seq_piece = ref_protein_seq_piece[0:len(ref_protein_seq_piece)-1]
                
                # make sure that the piece is location within the bounds of the exon translation    
                for a in list(re.finditer(str(ref_protein_seq_piece), str(ref_protein_seq))): 
                    if a.start() >= exon_start and a.end() <= exon_stop:
                        al_piece.set_protein_locations (a.start(), a.end())
                        break
       
            if al_piece.type == "insertion":
                al_piece.set_protein_locations(previous.ref_protein_stop, previous.ref_protein_stop + 1)
        
            previous = al_piece
            
    
    def determine_absolute_coordinates (self):
        '''
        Sets the absolute genomic locations for alignment pieces
        '''
            
        dmc = DataMapContainer.Instance ()
        conf_reader = ConfigurationReader.Instance ()
        
        expansion = int(conf_reader.get_value("local_ensembl", "expansion"))
        
        data_map = dmc.get((self.ref_protein_id, self.species))
        start = data_map.start
        alignment_start = self.alignment_exon.alignment_info["query_start"]
        
        for al_piece in self.alignment_pieces:
            
            if al_piece.type in ["coding", "insertion"]:
                
                real_start = max (1, start - expansion) + alignment_start + al_piece.alignment_start
                real_stop  = real_start + len(al_piece.ref_seq)
                
                al_piece.set_genomic_locations(real_start, real_stop, data_map.location_id)
                
                
    def create_cDNA (self):
        
        '''
        Create the cDNA for an alignment exon.
        Pad the gaps with Ns.
        '''
        total_exon_len = len(self.ref_exon.sequence)
        
        alignment_start = self.alignment_exon.alignment_info["sbjct_start"]
        padded_cdna = "N"* (alignment_start-1)
        
        len_added = len(padded_cdna)
        
        for al_piece in self.alignment_pieces:
            if al_piece.type == "coding":
                padded_cdna += al_piece.spec_seq
                len_added += len(al_piece.spec_seq)
            elif al_piece.type == "insertion":
                ns_to_add = (3 - (len(al_piece.spec_seq)) % 3) % 3
                padded_cdna += al_piece.spec_seq + "N"*ns_to_add
            elif al_piece.type == "deletion":
                padded_cdna += "N"*len(al_piece.spec_seq)
                len_added += len(al_piece.spec_seq)
                
        padded_cdna += "N" * (total_exon_len-len_added)
        return padded_cdna
                
            
                
        
        
                        
                