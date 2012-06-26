'''
Created on Apr 30, 2012

@author: ana
'''


# Python imports
import copy

# BioPython imports
from Bio                import SeqIO
from Bio.Seq            import Seq
from Bio.SeqRecord      import SeqRecord
from Bio.Alphabet       import IUPAC
from Bio.Alphabet.IUPAC import unambiguous_dna

# utilities imports
from utilities.DirectoryCrawler import DirectoryCrawler
from utilities.Logger           import Logger
from utilities.FileUtilities    import get_reference_species_dictionary,\
    read_seq_records_from_file

# data analysis imports
from data_analysis.utilities.ExonUtils          import LongestCommonSubstring
from data_analysis.base.EnsemblExon             import EnsemblExon
from data_analysis.containers.DataMapContainer  import DataMapContainer
from data_analysis.containers.ProteinContainer  import ProteinContainer
from data_analysis.containers.EnsemblExonContainer import EnsemblExonContainer

class EnsemblExons(object):
    '''
    Class which contains all the ensembl exons for a certain protein.
    '''

    def __init__(self, data_map_key, ref_species=None):
        '''
        @param data_map_key: (ref_protein_id, species)
        '''
        self.ref_protein_id    = data_map_key[0]
        self.species        = data_map_key[1]
        if not ref_species:
            spec_dict = get_reference_species_dictionary()
            ref_species = spec_dict[data_map_key[1]]
        self.ref_species    = ref_species
        self.exons          = {}
        
    def get_exon_file_path (self):
        '''
        Retrieve the file with the ensembl exons in fasta format
        '''
        dc = DirectoryCrawler()
        return "{0}/{1}.fa".format(dc.get_exon_ensembl_path(self.ref_protein_id), self.species)
    
    def load_exons (self):
        '''
        Load the exons from the fasta file and create
        the dictionary mapping them by their Ensembl id.
        Exons are given appropriate ordinals. 
        '''
        data_map_container      = DataMapContainer.Instance()
        logger                  = Logger.Instance()
        containers_logger       = logger.get_logger('containers')
        
        data_map = data_map_container.get((self.ref_protein_id, self.species))
        self.strand = data_map.strand
        
        fasta_path = self.get_exon_file_path()
        try:
            fasta = open(fasta_path, 'r')
        except IOError:
            containers_logger.error("%s,%s,%s" % (self.ref_protein_id, self.species, "Loading ensembl exons failed."))
            return None
         
        exon_list = []
        seq_records = read_seq_records_from_file(fasta, IUPAC.ambiguous_dna)
        
        for seq_record in seq_records:
            (start, stop, transcript_id, exon_id, strand) = seq_record.id.split('|')
            if (int(strand) == 1):
                self.strand = 1
                exon = EnsemblExon((self.ref_protein_id, self.species), exon_id, start, stop, strand, seq_record.seq)
            else:
                self.strand = -1
                exon = EnsemblExon((self.ref_protein_id, self.species), exon_id, stop, start, strand, seq_record.seq)
            exon_list.append(exon)
        fasta.close()
        self.exons = dict([(exon.exon_id, exon) for exon in exon_list])
        
        # assign orinals to exons
        ordinal = 1
        if self.strand == 1:
            for exon in sorted (self.exons.values(), key = lambda exon: exon.start ):
                exon.set_exon_ordinal(ordinal)
                ordinal += 1
        else:
            for exon in sorted (self.exons.values(), key = lambda exon: exon.start, reverse = True):
                exon.set_exon_ordinal(ordinal)
                ordinal += 1
        
        return exon_list
        
    def get_ordered_exons (self):
        '''
        Retrieves the exons in the correct protein coding order
        @return: list of ordered exons
        '''
        if self.strand == 1:
            return sorted (self.exons.values(), key = lambda exon: exon.start )
        else:
            return sorted (self.exons.values(), key = lambda exon: exon.start, reverse = True)
        
    def get_cDNA (self):
        '''
        @return: merged exon sequences
        '''
        exons = self.get_ordered_exons()
        merged_exons_seq = Seq("", IUPAC.ambiguous_dna)
        exon_locations = {}
        start = 1
        end = 1
        exon_id = 1
        for exon in exons:
            end += len(exon.sequence)-1
            exon_locations[exon_id] = (start, end)
            exon_id += 1
            start = end+1
            end = start
            merged_exons_seq += exon.sequence
        cdna_seq = SeqRecord(seq=merged_exons_seq, id=self.species, description="")
        return (cdna_seq, exon_locations)
    
    def get_coding_cDNA(self):
        '''
        @return: the cDNA without the UTR regions
        '''
        dmc = DataMapContainer.Instance()
        pc = ProteinContainer.Instance()
        data_map = dmc.get((self.ref_protein_id, self.species))
        
        protein = pc.get(data_map.protein_id)
        protein_sequence = str(protein.get_sequence_record().seq)
        
        (cDNA_record, locs) = self.get_cDNA()
        cDNA = cDNA_record.seq
        longest_frame = -1
        longest_translation = ""
        for frame in (0,1,2):
            translation = cDNA[frame:].translate()
            longest_substring = LongestCommonSubstring (translation, protein_sequence)
            if len(longest_substring) > len(longest_translation):
                longest_translation = longest_substring
                actual_translation = translation
                longest_frame = frame
            
        prot_start = actual_translation.find(longest_translation)
        cDNA_start = prot_start*3 + longest_frame
        cDNA_end = (prot_start + len(longest_translation) + 1) * 3 + longest_frame
        
        return (cDNA[cDNA_start:cDNA_end], cDNA_start, cDNA_end)

    
    def export_coding_exons_to_fasta (self, fasta_file):
        '''
        Exports the protein coding exons to .fasta file
        @param fasta_file: the output fasta file
        '''
        dmc = DataMapContainer.Instance()
        data_map = dmc.get((self.ref_protein_id, self.species))
        
        new_exons = self.get_coding_exons()
        exon_records = []
        # >969067|969174|ENSAMET00000013141|ENSAMEE00000125733|1
        for exon in new_exons:
            exon_id = "%d|%d|%s|%s|%d" % (exon.start, exon.stop, data_map.transcript_id, exon.exon_id, exon.strand)
            record = SeqRecord(seq = exon.sequence, id = exon_id, description = "")
            exon_records.append(record)
            
        SeqIO.write(exon_records, fasta_file, "fasta")
        
    def get_coding_exons (self):
        
        (coding_cDNA, cDNA_start, cDNA_end) = self.get_coding_cDNA()
        
        total_start = 0
        total_end = 0
        ordered_exons = self.get_ordered_exons()
        coding_exons = []
        
        for exon in self.get_ordered_exons():
            
            total_start = total_end 
            total_end += len(exon.sequence)

            if total_end <= cDNA_start:
                continue
            if total_start >= cDNA_end:
                continue
            
            if total_start >= cDNA_start and total_end <= cDNA_end:
                coding_exon = copy.deepcopy(exon)
                
            if total_start < cDNA_start and total_end > cDNA_start and total_end <= cDNA_end:
                new_start = cDNA_start - total_start
                coding_exon = copy.deepcopy(exon)
                coding_exon.sequence = coding_exon.sequence[new_start:]
                coding_exon.start += new_start
                
            if total_start < cDNA_end and total_end > cDNA_end and total_start >= cDNA_start:
                coding_exon = copy.deepcopy(exon)           
                new_end = len(coding_exon.sequence) - (total_end - cDNA_end)
                coding_exon.end = new_end
                coding_exon.sequence = coding_exon.sequence[:new_end]
                coding_exon.end -= (total_end - cDNA_end)
                
            if total_start < cDNA_start and total_end > cDNA_end:
                coding_exon = copy.deepcopy(exon)
                new_start = cDNA_start - total_start
                new_end = len(coding_exon.sequence)- (total_end - cDNA_end)
                coding_exon.sequence = coding_exon.sequence[new_start:new_end]
                coding_exon.start += new_start
                coding_exon.stop -= (total_end - cDNA_end)
                
                
            coding_exons.append(coding_exon)
            
        return coding_exons
    
    def get_exon_ids_from_ccDNA_locations (self, start, stop):
        '''
        Retrieves the exons (with their relative locations set)
        for certain locations on the coding cDNA.
        '''
        # minor correction due to the unnatural
        # way of marking the beginning the alignment outputs
        start -= 1
        
        coding_exons = self.get_coding_exons()
        #coding_exons = self.get_ordered_exons()
        start_loc_on_genome = coding_exons[0].start
        present_exons = []
        relative_start,relative_stop = 0,0
        
        for exon in coding_exons:
            # refresh the start to previous stop
            relative_start = relative_stop
            relative_stop += len(exon.sequence)
            #(relative_start,relative_stop) = (exon.start - start_loc_on_genome, exon.stop - start_loc_on_genome)
            # ------|----------|--------
            # --------------------------
            if relative_stop <= start or relative_start >= stop:
                continue
            # ------|----------|--------
            #---------|---|-------------
            if relative_start >= start and relative_stop <= stop:
                exon.set_relative_start (0)
                exon.set_relative_stop (relative_stop - relative_start)
                present_exons.append(exon)
                continue
            # ------|----------|--------
            # ---|----------------|-----
            if relative_start <= start and relative_stop >= stop:
                exon.set_relative_start (start - relative_start)
                exon.set_relative_stop (stop - relative_start)
                present_exons.append(exon)
                continue
            
            # ------|----------|--------
            # ---|--------|-------------
            if relative_start < start and relative_stop <= stop:
                exon.set_relative_start (start - relative_start)
                exon.set_relative_stop (relative_stop - relative_start)
                present_exons.append(exon)
                
            # ------|----------|--------
            # --------|-----------|-----
            if relative_start >= start and relative_stop > stop:
                exon.set_relative_start(0)
                exon.set_relative_stop (stop - relative_start)
                present_exons.append(exon)
        
        return present_exons
    
    
    def set_coding_exon_frames (self):
        '''
        Serves to set the frames for coding exons.
        The frame of the first coding exon is always 0. The others
        are determined based on the sum of the previous exons lengths. 
        '''
        
        coding_exons = self.get_coding_exons()
        ec = EnsemblExonContainer.Instance()
        length = 0
        
        for exon in coding_exons:
            ensembl_exon = ec.get(exon.exon_id)
            if length % 3 == 0:
                ensembl_exon.set_frame (0) 
                self.exons[exon.exon_id].set_frame (0) 
                exon.set_frame (0)
            elif length % 3 == 1:
                ensembl_exon.set_frame (2)
                self.exons[exon.exon_id].set_frame (2) 
                exon.set_frame (2)
            else:
                ensembl_exon.set_frame (1)
                self.exons[exon.exon_id].set_frame (1)
                exon.set_frame (1) 
                
            length += len(exon.sequence)
            
        return coding_exons
        
# ana i luka zaljubljeni par
# kad se budu zenili dobiceju dar    

                
 
if __name__ == '__main__':
    pass
    
    
        
        