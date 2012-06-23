'''
Created on Jun 21, 2012

@author: intern
'''
from data_analysis.containers.ExonContainer import ExonContainer
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
from data_analysis.containers.EnsemblExonContainer import EnsemblExonContainer
from data_analysis.containers.ProteinContainer import ProteinContainer

def translate_ensembl_exons (ensembl_exons):
    
    cdna = ""
    #determine the frame
    first_exon = ensembl_exons[0]
    if first_exon.relative_start == 0:
        new_frame = first_exon.frame
    else:
        frame_complement = (first_exon.relative_start - first_exon.frame) % 3
        if frame_complement == 1:
            new_frame = 2
        elif frame_complement == 2:
            new_frame = 1
        else:
            new_frame = 0 
        
    #create the cdna
    for ee in ensembl_exons:
        new_sequence = ee.sequence[ee.relative_start:ee.relative_stop]
        cdna += new_sequence
    al_exon = AlignmentExonPiece("ensembl", -1, "", cdna)
    al_exon.set_frame(new_frame)
    return cdna[new_frame:].translate()


def process_insertion_in_exon(human_seq, spec_seq, ordinal):
    '''
    Annotate insertion. Insertion can be a 'real' insertion, or 
    marked as a frameshift (meaning that it is 1, 2, 4 or 5 
    base pairs long)
    @return: AlignmentExonPiece
    '''
    if len(spec_seq) in (1,2,4,5):
        exon_type = "frameshift"
    else:
        exon_type = "insertion"
        
    al_exon = AlignmentExonPiece(exon_type, ordinal, human_seq, spec_seq)
    return al_exon
        


def process_insertion_free_region(human_seq, spec_seq, frame, ordinal):
    '''
    Process a region of an alignment without any insertions in the 
    reference species exon. There can exist deletions however. 
    '''
    al_exons = []
    
    # determine the number of coding sequences interspersed with gaps 
    number_of_ungapped_sequences = len(re.split("-+", str(spec_seq)))
    # create a pattern to load the sequences
    pattern_string = "([ATGCN]+)"
    for i in range (0, number_of_ungapped_sequences-1):
        pattern_string += "([-]+)([ATGCN]+)"
    pattern = re.compile (pattern_string)
    sequences = re.match (pattern, str(spec_seq))
    
    
    # auxiliary variables
    in_coding = True
    start, stop = 0,0
    coding_len = 0
    
    for seq in sequences.groups():
        stop += len(seq)
        
        species = seq
        human   = human_seq[start:stop]
        
        # if we're not in the coding region, there is nothing to process
        # non-coding here means deletion in the referent exon
        if in_coding:
            if coding_len == 0:
                new_frame = frame
            else:
                frame_status = abs(coding_len - frame) % 3
                if frame_status == 1:
                    new_frame = 2
                elif frame_status == 2:
                    new_frame = 1
                else:
                    new_frame = 0
                
            exon = AlignmentExonPiece("coding", ordinal, human, species)
            exon.set_frame(new_frame)
            al_exons.append(exon)
            
        # there is a deletion. Remember this also - if the deletion is
        # one of two bases long, we pad it with Ns and translate regularly    
        else:
            exon = AlignmentExonPiece("deletion", ordinal, human, species)
            al_exons.append(exon)
            
            
        in_coding = not in_coding
        start = stop
        coding_len += len(seq)
            
    return al_exons


def split_exon_seq (alignment_exon, coding_exon):
    '''
    Splits  the exon alignment into pieces: the insertions in the
    referent exon and the parts without insertions (but with
    possible deletions)
    @return: all the alignment exon pieces for this exon alignment
    '''
    
    genomic_dna = alignment_exon.alignment_info["query_seq"]
    exon_dna    = alignment_exon.alignment_info["sbjct_seq"]
    
    number_of_ungapped_sequences = len(re.split("-+", str(exon_dna)))
    pattern_string = "([ATGCN]+)"
    for i in range (0, number_of_ungapped_sequences-1):
        pattern_string += "([-]+)([ATGCN]+)"
    pattern = re.compile (pattern_string)
    sequences = re.match (pattern, str(exon_dna))
    
    alignment_pieces = []
    

    # determine starting frame (in need of refactoring)
    if alignment_exon.alignment_info["sbjct_start"] == 1:
        frame = coding_exon.frame
    elif alignment_exon.alignment_info["sbjct_start"] == 2:
        if coding_exon.frame == 0:
            frame = 2
        elif coding_exon.frame == 1:
            frame = 0
        else:
            frame = 1
    else:
        frame_status = abs(alignment_exon.alignment_info["sbjct_start"] - 1 - coding_exon.frame) % 3
        if frame_status == 2:
            frame = 1
        elif frame_status == 1:
            frame = 2
        else:
            frame = 0
    
    ## auxiliary variables ##        
    coding_len = 0
    start, stop = 0,0
    ordinal = 0
    in_coding = True
            
    ########## PROCESS SEQUENCES (INSERTIONS / CODING REGIONS) ############
    
    for seq in sequences.groups ():
        
        stop += len(seq)
        
        # determine the working sequences
        human_seq = seq
        spec_seq  = genomic_dna[start:stop]
        
        ## if not coding, then we have an insertions ##
        if not in_coding:
            exon = process_insertion_in_exon (human_seq, spec_seq, ordinal)
            alignment_pieces.append(exon)

        # determine frame
        if in_coding:
            if coding_len == 0:
                new_frame = frame
            else:
                frame_status = abs(coding_len - frame) % 3
                if frame_status == 1:
                    new_frame = 2
                elif frame_status == 2:
                    new_frame = 1
                else:
                    new_frame = 0
        
        # if we're in the coding region, process it further for deletions    
        if in_coding:
            al_exons = process_insertion_free_region (human_seq, spec_seq, new_frame, ordinal)
            alignment_pieces.extend(al_exons)
        
        
        # update the status' and lengths
        start = stop
        if in_coding:
            coding_len += len(seq)
        in_coding = not in_coding
        ordinal += 1
        
    return alignment_pieces
        
    
    
def patch_sequences(alignment_exon, actual_exon_len):  
    '''  
    Pads the alignment sequences with leftovers from the previous    
    exon, with the right amount of Ns in the beginning (for the initial gaps)
    and the right amount of Ns in the end of the sequence   
    '''     
    genomic_cdna = Seq("", IUPAC.ambiguous_dna)  
    exon_cdna    = Seq("", IUPAC.ambiguous_dna)      
    # decide what to append to the cDNA   
    (translation_start, translation_end) = (alignment_exon.alignment_info["sbjct_start"],   
                                            alignment_exon.alignment_info["sbjct_end"])
  
    genomic_cdna +=  ( "N"* (translation_start - 1)) 
    exon_cdna    +=  ("N" * (translation_start - 1)) 
 
    genomic_cdna += alignment_exon.alignment_info["query_seq"]  
    exon_cdna    += alignment_exon.alignment_info["sbjct_seq"]
    
    genomic_cdna += "N" * (actual_exon_len - translation_end)   
    exon_cdna    += "N" * (actual_exon_len - translation_end)    
 
    return (genomic_cdna, exon_cdna)   
    



def translate_alignment_exons(ref_protein_id, reference_species, alignment_exon):
    
    # load the reference exons
    exon_container  = ExonContainer.Instance()
    reference_exons = exon_container.get ((ref_protein_id, reference_species, "ensembl"))
    
    # get just the coding exons
    for coding_exon in reference_exons.get_coding_exons():
        if coding_exon.exon_id == alignment_exon.ref_exon_id:
            break
    coding_seq_len = len(coding_exon.sequence)
    
    
    (genomic_cdna, exon_cdna) = patch_sequences (alignment_exon, coding_seq_len)
    
    return (genomic_cdna, exon_cdna)
    
    
class BestExonAlignment (object):
    '''
    Class containing all the alignments for one exon, with a status
    describing which of the alignments is the best
    '''
    def __init2__ (self, ref_exon_id, gene_alignment_exon, cdna_alignment_exon, ensembl_exons):
        pass
    def __init__ (self, ref_exon_id, sw_gene_exon = None, ensembl_exons = None, status = None, ensembl_alignment = None, sw_gene_alignment = None):
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
        print "whole exon:", self.ref_exon.exon_id, self.ref_exon.frame , complete_protein_exon_seq
        exon_prot_start = str(ref_protein_seq).find(str(complete_protein_exon_seq))
        exon_prot_stop = exon_prot_start + len(complete_protein_exon_seq)
        
        for frame in (0,1,2):
            partial_protein_ref_seq = partial_ref_seq[frame:].translate()
            if str(partial_protein_ref_seq).endswith("*"):
                partial_protein_ref_seq = partial_protein_ref_seq[0:len(partial_protein_ref_seq)-1]
            print partial_protein_ref_seq
            found = False
            if str(complete_protein_exon_seq).find (str(partial_protein_ref_seq)) != -1:
                for a in list(re.finditer(str(partial_protein_ref_seq), str(ref_protein_seq))): 
                    print a.start(), a.end()
                    if a.start() >= exon_prot_start and a.end() <= exon_prot_stop:
                        self.ref_protein_seq = partial_protein_ref_seq
                        self.ref_protein_start = a.start()
                        self.ref_protein_stop = a.end()
                        found = True
                        break
            if found:
                break
                    
        
        self.spec_protein_seq = translate_ensembl_exons(self.ensembl_exons)
        
        
   
class SWGeneAlignment (object):
    
    def __init__ (self, ref_protein_id, ref_exon, alignment_exon):   
        self.ref_protein_id = ref_protein_id
        self.ref_exon = ref_exon
        self.alignment_exon = alignment_exon
        self.load_alignment_pieces()  
        
    def load_alignment_pieces (self):
        self.alignment_pieces = split_exon_seq(self.alignment_exon, self.ref_exon)
        

class AlignmentExonPiece (object):
    '''
    Class representing one alignment piece. It is a part 
    of one exon alignment. It can be of various types:
     - coding: means that it corresponds to the coding part of the reference species DNA
     - insertion: meaning that it is an insertion in the coding ref species DNA
     - frameshift: meaning that is an insertion, but of length 1,2,4 or 5 base pairs
    '''
    def __init__ (self, exon_type, ordinal, ref_seq, spec_seq):
        '''
        @param exon_type: coding / insertion / frameshift
        @param ref_seq: nucleotide sequence in the reference species
        @param spec_seq: nucleotide sequence in the species
        @param ordinal: ordinal in the one exon alignment
        '''
        self.type       = exon_type
        self.ordinal    = ordinal
        self.ref_seq    = ref_seq
        self.spec_seq   = spec_seq
        
    def set_location (self, start, stop):
        '''
        Set location relative to the alignment
        '''
        self.start = start
        self.stop = stop
        
    def set_frame(self, frame):
        '''
        Set the translation frame (0,1,2)
        '''
        self.frame = frame
        
    def set_translations (self, ref_spec_protein_seq, spec_protein_seq):
        '''
        Set the protein translation of the alignment piece.
        These alignments need not completely correspond to the DNA sequence, 
        because some of them might be extended on both sides to account 
        for the translation frame (so as not to lose AAs)
        '''
        self.ref_spec_protein_seq = ref_spec_protein_seq
        self.spec_protein_seq = spec_protein_seq
        
    def set_protein_locations (self, ref_protein_start, ref_protein_stop):
        '''
        Set the locations of the translated piece on the 
        referent species protein
        '''
        self.ref_protein_start = ref_protein_start
        self.ref_protein_stop = ref_protein_stop
        
    
        








