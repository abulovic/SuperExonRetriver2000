'''
Created on Mar 22, 2012

@author: marioot
'''
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re, sys
from operator import attrgetter
from timeit import itertools
from utilities.DescriptionParser import DescriptionParser
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
from data_analysis.containers.ProteinContainer import ProteinContainer

def get_exon(exons, exon_id):
    for exon in exons:
        if exon.id is exon_id:
            return exon
        
def set_exon(exons, exon):
    for exon_id in exons:
        if exon.id is exon_id:
            exon_id = exon
    return exons

def find_frame(coding_dna, target_transcript):
    dna = Seq(coding_dna, IUPAC.unambiguous_dna)
    for frame in range(3):
        translation             = dna[frame:].translate().tostring()
        location_on_transcript  = target_transcript.find(translation)
        if location_on_transcript is not -1:
            return [frame, 
                    translation, 
                    target_transcript[location_on_transcript + len(translation):]]

def indentify_initial_exon_frames(exons, target_prot):
    for exon in exons:
        frame = find_frame(exon.target.replace("-", ""), 
                           target_prot)[0]
        exon.set_frame(frame)
    return exons

def repair_broken_exons(exons, viable_exons):
    for i in range(len(viable_exons) - 1):
        if viable_exons[i + 1] == viable_exons[i] + 1:
            exon_ = get_exon(exons, viable_exons[i])
            _exon = get_exon(exons, viable_exons[i + 1])
            if exon_.t_end == exon_.length and _exon.t_start == 1 and _exon.frame != 0:
                query_nucleotides   =   exon_.query[len(exon_.query) - (3 - _exon.frame) : ]
                target_nucleotides  =   exon_.target[len(exon_.target) - (3 - _exon.frame) : ]
                
                _exon.query         =   query_nucleotides + _exon.query
                _exon.target        =   target_nucleotides + _exon.target
                
                exon_.query         =   exon_.query[ : len(exon_.query) - (3 - _exon.frame)]
                exon_.target        =   exon_.target[ : len(exon_.target) - (3 - _exon.frame)]
                
                exons               =   set_exon(exons, exon_)
                exons               =   set_exon(exons, _exon)
           
def translate_dna(coding_dna, frame):
    dna = Seq(coding_dna, IUPAC.unambiguous_dna)
    return dna[frame:].translate().tostring()



def get_to_next_valid_position(dna, position):
    while dna[position] is '-':
        position    += 1
        if position == len(dna): 
            return -1
    return position

def round_up_to_triplet(position):
    if (position + 1) % 3 == 0:
        return position + 1
    else:
        return position + 2

def count_initial_gaps(target_dna):
    cnt = 0
    while target_dna[cnt] is '-':
        cnt += 1
    return cnt
    
def analyse_query_dna(query_dna, target_dna, target_prot):
    query_prot  =   ""
    init_gaps   =   count_initial_gaps(query_dna)
    query_dna   =   query_dna[init_gaps:]
    target_dna  =   target_dna[init_gaps:]
    
    next_gap    =   query_dna.find('-')
    if next_gap is not -1:
        [frame, left_target_prot, right_target_prot] = find_frame(target_dna[:next_gap], 
                                                                  target_prot) 
        if frame > 0:
            query_prot += "X"
        query_prot          +=  translate_dna(query_dna[:next_gap], 
                                              frame)
        next_valid_position =   get_to_next_valid_position(query_dna, 
                                                           next_gap)
        number_of_gaps      =   next_valid_position - next_gap
        untranslated_tail   =   len(query_dna[frame:next_gap])%3
        if (untranslated_tail > 0) and (untranslated_tail + number_of_gaps) % 3 is 0:
            query_prot      += "X"
            
        query_prot          +=  analyse_query_dna(query_dna[next_valid_position:], 
                                                  target_dna[next_valid_position:], 
                                                  right_target_prot)
        query_dna           =   query_dna[:next_gap]
    else:
        frame               =   find_frame(target_dna, 
                                           target_prot)[0]
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  translate_dna(query_dna, 
                                              frame)
    return query_prot
    
def analyse_SW_alignment(query_dna, target_dna, target_prot):
    query_prot  =   ""
    init_gaps   =   count_initial_gaps(target_dna)
    insertion   =   translate_dna(query_dna[0:init_gaps], 0)
    query_prot  +=  insertion
    query_dna   =   query_dna[init_gaps:]
    target_dna  =   target_dna[init_gaps:]
    next_gap    =   target_dna.find('-')
    if next_gap is not -1:
        #not the end of sequence
        [frame, left_target_prot, right_target_prot] = find_frame(target_dna[:next_gap], 
                                                                  target_prot)        
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  analyse_query_dna(query_dna[frame:next_gap], 
                                                  target_dna[frame:next_gap], 
                                                  left_target_prot)
        next_valid_position =   get_to_next_valid_position(target_dna, 
                                                           next_gap)
        number_of_gaps      =   next_valid_position - next_gap
        insertion_size      =   number_of_gaps/3
        if insertion_size > 0:
            query_prot      +=  translate_dna(query_dna[next_gap:next_gap + 3*insertion_size], 0).lower()
        query_prot          +=  analyse_SW_alignment(query_dna[next_valid_position:], 
                                                     target_dna[next_valid_position:], 
                                                     right_target_prot)
    else:
        #end of the sequence
        [frame, left_target_prot, target_prot] = find_frame(target_dna, 
                                                            target_prot)        
        if frame > 0:
            query_prot      +=  "X"
        query_prot          +=  analyse_query_dna(query_dna[frame:], 
                                                  target_dna[frame:], 
                                                  left_target_prot)
        #check to see if the last segment ends as a uncomplete triplet
        if len(target_dna[frame:])%3 > 0:
            query_prot      +=  "X"
    return query_prot

def transcribe_exons(exons, target_prot):
    query_prot      = ""
    viable_exons    = []
    for exon in exons:
        if exon.viability:
            viable_exons.append(exon.id)
    indentify_initial_exon_frames(exons, 
                                  target_prot)
    viable_exons = sorted(viable_exons)
    repair_broken_exons(exons, 
                        viable_exons)
    for exon in exons:
        if exon.viability:
            query_prot += analyse_SW_alignment(exon.query, 
                                               exon.target, 
                                               target_prot)
    return query_prot

def transcribe_exons2 (protein_id, species, alignment_alg):
    
    pc = ProteinContainer.Instance()
    protein = pc.get(protein_id)
    
    print protein.get_sequence_record()
    

def transcribe_aligned_exons (protein_list, algorithms):
    
    for protein_id in protein_list:
        # get list of available species:
        species_list = DescriptionParser().get_species(protein_id)
        for species in species_list:
            for alignment_alg in algorithms:
                
                transcribe_exons2 (protein_id, species, alignment_alg)