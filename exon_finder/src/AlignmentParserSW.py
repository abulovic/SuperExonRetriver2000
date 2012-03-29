'''
Created on Mar 27, 2012

@author: Mario, anana
'''
from AlignmentParser import AlignmentParser
import re, os
from Exon import Exon
from operator import attrgetter
from timeit import itertools

class AlignmentParserSW (AlignmentParser):



    def __init__(self, algorithmName="Smith-Waterman"):
        '''
        Constructor
        '''
        self.algorithmName = algorithmName
        super(AlignmentParserSW, self).__init__()
        
    def batchParseOutputExDna (self):
        
        swoutDirectory = "%s/dna/" % (self.swoutAbs)
        exonsDir = {}
        
        for swoutFile in os.listdir(swoutDirectory):
            if (swoutFile.endswith('.swout')):
                species = swoutFile.split('.')[0]
                exons = self.parseOutputSW(swoutFile)
                exonsDir[species] = exons
                
        return exonsDir
                
        
    def parseOutputSW(self, swout_path):
        exons = []
        swout = open(swout_path, 'r')
        exon_block_pattern  = re.compile(r'>(\d+) exon length (\d+)')
        intervals_pattern   = re.compile(r'Intervals: (\d+) (\d+) (\d+) (\d+)')
        identity_pattern    = re.compile(r'Identity: (\d+)/(\d+) \(\d+\.\d+%\)')
        score_pattern       = re.compile(r'Score: (\d+)\.\d+')
        sequence_pattern    = re.compile(r'[\d\s]+([ATCG-]+)[\d\s]+')
        query_sequence_flag = True
        for line in swout.readlines():
            exon_block_match    =   re.match(exon_block_pattern, line)
            if exon_block_match is not None:
                exons.append(Exon(exon_block_match.groups()[0], 
                                  exon_block_match.groups()[1], 
                                  "", 
                                  ""))
                continue
            intervals_match     =   re.match(intervals_pattern, line)
            if intervals_match  is not None:
                exons[len(exons) - 1].set_intervals(intervals_match.groups()[0], 
                                                    intervals_match.groups()[1], 
                                                    intervals_match.groups()[2], 
                                                    intervals_match.groups()[3])
                continue
            identity_match      =   re.match(identity_pattern, line)
            if identity_match   is not None:
                exons[len(exons) - 1].set_identity(identity_match.groups()[0],
                                                   identity_match.groups()[1])
                continue
            score_match         =   re.match(score_pattern, line)
            if score_match      is not None:
                exons[len(exons) - 1].set_score(score_match.groups()[0])
                continue
            sequence_match      =   re.match(sequence_pattern, line)
            if sequence_match   is not None:
                if query_sequence_flag is True:
                    exons[len(exons) - 1].query     +=  sequence_match.groups()[0]
                else:
                    exons[len(exons) - 1].target    +=  sequence_match.groups()[0]
                query_sequence_flag                 =   not query_sequence_flag
                continue
        return exons
    
    def isSorted(self, inList):
        for i in range(len(inList) - 1):
            if inList[i] > inList[i+1]: return False
        return True     
    
    def calculate_total_score(self, valid_id_list, exons):
        score = 0.0
        for exon in exons:
            if exon.id in valid_id_list:
                score += exon.fitness()
        return score
    
    def exaustive_search(self, exons):
        full_array          = []
        highest_score       = 0.0
        best_combination    = []
        for exon in exons:
            full_array.append(exon.id)
        for i in range(len(full_array)):  
            all_combinations                = itertools.combinations(full_array, len(full_array) - i)
            for comb in all_combinations:
                if self.isSorted(comb): 
                    score                   = self.calculate_total_score(comb, exons)
                    if score > highest_score: 
                        highest_score       = score
                        best_combination    = comb
            if highest_score > 0.0: return best_combination
    
    def discard_FP(self, exons):
        sorted_exons    = sorted(exons, key = attrgetter('q_start', 'q_end'))
        best_solution   = self.exaustive_search(sorted_exons)
        for exon in exons:
            if exon.id in best_solution:
                exon.set_viablity(True)
            else:
                exon.set_viablity(False)
        return exons

    
if __name__ == '__main__':
    swParser = AlignmentParserSW()
    swParser.setProteinFolder("ENSP00000311134")
    exons = swParser.parseOutput("/home/intern/Documents/sw_input/sw_output/sus_scrofa_out.txt")
    for exon in sorted(exons, reverse=True):
        print exon.id, float(exon.no_of_matches) / exon.alignment_length

        