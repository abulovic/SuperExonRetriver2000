'''
Created on Jun 21, 2012

@author: anana

Auxiliary functions to utilize other classes in dna to protein translation.
'''
# Python imports
import re

# BioPython imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# data analysis imports
from data_analysis.translation.AlignmentExonPiece import AlignmentExonPiece

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


def process_insertion_in_exon(human_seq, spec_seq, ordinal, frame):
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
    al_exon.set_frame(frame)
    return al_exon
        
def replace_slash (matchobj):
    (ch1, slash, ch2) = matchobj.groups()
    if slash == "-":
        return "%sN%s" % (ch1, ch2)
    elif slash == "--":
        return "%sNN%s" % (ch1, ch2)

def process_insertion_free_region(human_seq, spec_seq, frame, ordinal, alignment_start):
    '''
    Process a region of an alignment without any insertions in the 
    reference species exon. There can exist deletions however. 
    If the deletions are 1 or 2 bases long, they are replaced
    by N (regarded as an assembly error)
    '''
    al_exons = []
    
    spec_seq = re.sub("([ATGCN])(-{1,2})([ATGCN])", replace_slash, spec_seq)
    
    # determine the number of coding sequences interspersed with gaps 
    number_of_ungapped_sequences = len(re.split("-+", str(spec_seq)))
    if number_of_ungapped_sequences == 1:
        al_piece = AlignmentExonPiece("coding", ordinal, human_seq, spec_seq)
        al_piece.set_frame(frame)
        al_piece.set_alignment_locations(alignment_start, alignment_start + len(spec_seq))
        return [al_piece]
    
    # create a pattern to load the sequences
    pattern_string = "([ATGCN]+)"
    for i in range (0, number_of_ungapped_sequences-1):
        pattern_string += "(-+)([ATGCN]+)"
    pattern = re.compile (pattern_string)
    sequences = re.match (pattern, str(spec_seq))
    
    
    # auxiliary variables
    in_coding = True
    start, stop = 0,0
    coding_len = 0
    alignment_stop = alignment_start
    
    for seq in sequences.groups():
        stop += len(seq)
        alignment_stop += len(seq)
        
        species = seq
        human   = human_seq[start:start + len(seq)]
        
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
            exon.set_alignment_locations(alignment_start, alignment_stop)
            al_exons.append(exon)
            
        # there is a deletion. Remember this also - if the deletion is
        # one of two bases long, we pad it with Ns and translate regularly    
        else:
            exon = AlignmentExonPiece("deletion", ordinal, human, species)
            al_exons.append(exon)
            
            
        in_coding = not in_coding
        start = stop
        alignment_start = alignment_stop
        coding_len += len(seq)
        ordinal += 1
            
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
        
        ## if not coding, then we have an insertions ##
        if not in_coding:
            exon = process_insertion_in_exon (human_seq, spec_seq, ordinal, new_frame)
            exon.set_alignment_locations(start, stop)
            alignment_pieces.append(exon)
        
        # if we're in the coding region, process it further for deletions    
        if in_coding:
            al_exons = process_insertion_free_region (human_seq, spec_seq, new_frame, ordinal, start)
            ordinal += len(al_exons)
            alignment_pieces.extend(al_exons)
        
        
        # update the status' and lengths
        start = stop
        if in_coding:
            coding_len += len(seq)
        in_coding = not in_coding
        ordinal += 1
        
    return alignment_pieces


    

def set_protein_sequences (alignment_pieces):
    '''
    Set the protein sequences for alignment pieces.
    For each of the coding alignment pieces, a check is performed:
    if there are coding pieces to the left and the right, the protein
    sequences are expanded to allow for the missing AAs due to 
    'middle of the codon' alignment start or stop
    '''
    
    for i in range (0, len(alignment_pieces)):
        
        al_piece = alignment_pieces[i]

        if al_piece.type == "insertion":
            
            add_beg_ref_seq, add_beg_spec_seq = "", ""
            add_end_ref_seq, add_end_spec_seq = "", ""
            
            # if there exists a left alignment piece
            if i-1 in range (0, len(alignment_pieces)):
                prev_piece          = alignment_pieces[i-1]
                how_much_to_take    = (3 - al_piece.frame) % 3
                add_beg_spec_seq    = prev_piece.spec_seq[len(prev_piece.spec_seq)-how_much_to_take:len(prev_piece.spec_seq)]
                
            # if there exists a right alignment piece    
            if i+1 in range (0, len(alignment_pieces)):
                next_piece          = alignment_pieces[i+1]
                how_much_to_take    = (3 - (len(al_piece.ref_seq) - al_piece.frame)) % 3
                add_end_spec_seq    = next_piece.spec_seq[0:how_much_to_take]

            spec_seq_to_translate = Seq (add_beg_spec_seq + al_piece.spec_seq + add_end_spec_seq, IUPAC.ambiguous_dna)
            
            if not add_beg_ref_seq:
                spec_protein_seq = spec_seq_to_translate[al_piece.frame:].translate()
            else:
                spec_protein_seq = spec_seq_to_translate.translate()
                
            al_piece.set_translations (None, spec_protein_seq)
            
                    
        if al_piece.type == "coding":
            
            add_beg_ref_seq, add_beg_spec_seq = "", ""
            add_end_ref_seq, add_end_spec_seq = "", ""
            
            # if there is a left alignment piece
            if i-1 in range (0, len(alignment_pieces)):
                
                prev_piece           = alignment_pieces[i-1]
                if prev_piece.type in ["frameshift", "insertion"]:
                    prev_piece       = alignment_pieces[i-2]
                    how_much_to_take = (3 - al_piece.frame) % 3
                    add_beg_ref_seq  = prev_piece.ref_seq [len(prev_piece.ref_seq)-how_much_to_take:len(prev_piece.ref_seq)]
                    add_beg_spec_seq = prev_piece.spec_seq[len(prev_piece.spec_seq)-how_much_to_take:len(prev_piece.spec_seq)]
                    
            # if there is a right alignment piece
            if i+1 in range (0, len(alignment_pieces)):
                
                next_piece           = alignment_pieces[i+1]
                if next_piece.type in ["frameshift", "insertion"]:
                    next_piece       = alignment_pieces[i+2]
                    how_much_to_take = (3 - (len(al_piece.ref_seq) - al_piece.frame)) % 3
                    add_end_ref_seq  = next_piece.ref_seq [0:how_much_to_take]
                    add_end_spec_seq = next_piece.spec_seq[0:how_much_to_take]
                    
            ref_seq_to_translate     = Seq (add_beg_ref_seq + al_piece.ref_seq + add_end_ref_seq, IUPAC.ambiguous_dna)
            spec_seq_to_translate    = Seq (add_beg_spec_seq + al_piece.spec_seq + add_end_spec_seq, IUPAC.ambiguous_dna)
            
            if not add_beg_ref_seq:
                ref_protein_seq      = ref_seq_to_translate[al_piece.frame:].translate()
                spec_protein_seq     = spec_seq_to_translate[al_piece.frame:].translate()
            else:
                ref_protein_seq      = ref_seq_to_translate.translate()
                spec_protein_seq     = spec_seq_to_translate.translate()
            
            al_piece.set_translations (ref_protein_seq, spec_protein_seq)
            
    return alignment_pieces
            
            
        


        
    
        








