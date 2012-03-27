'''
Created on Mar 27, 2012

@author: Mario, anana
'''
from AlignmentParser import AlignmentParser
import re
from Exon import Exon

class swAlignmentParser(AlignmentParser):



    def __init__(self, algorithmName="Smith-Waterman"):
        '''
        Constructor
        '''
        self.algorithmName = algorithmName
        super(swAlignmentParser, self).__init__()
        
    def parseOutput (self, alignmentOutputFile, numberOfExons):
        '''
        @return: list of exons 
        '''
        exons = []
        swout = open(alignmentOutputFile, 'r')
        exon_block_pattern = re.compile(r'>(\d+) exon length (\d+)')
        intervals_pattern =  re.compile(r'Intervals: (\d+) (\d+) (\d+) (\d+)')
        identity_pattern =   re.compile(r'Identity: (\d+)/(\d+) \(\d+\.\d+%\)')
        score_pattern =      re.compile(r'Score: (\d+)\.\d+')
        sequence_pattern =   re.compile(r'[\d\s]+([ATCG-]+)[\d\s]+')
        
        query_sequence_flag = True
        for line in swout.readlines():
            
            exon_block_match =  re.match(exon_block_pattern, line)
            if exon_block_match is not None:
                exons.append(Exon(exon_block_match.groups()[0], "", "", int(exon_block_match.groups()[1])))
                continue
            
            intervals_match =   re.match(intervals_pattern, line)
            if intervals_match  is not None:
                exons[len(exons) - 1].set_intervals(intervals_match.groups()[0], 
                                                    intervals_match.groups()[1], 
                                                    intervals_match.groups()[2], 
                                                    intervals_match.groups()[3])
                continue
            
            identity_match =    re.match(identity_pattern, line)
            if identity_match   is not None:
                exons[len(exons) - 1].set_identity(identity_match.groups()[0],
                                                   identity_match.groups()[1])
                continue
            
            score_match =       re.match(score_pattern, line)
            if score_match      is not None:
                exons[len(exons) - 1].set_score(score_match.groups()[0])
                continue
            
            sequence_match =    re.match(sequence_pattern, line)
            if sequence_match   is not None:
                if query_sequence_flag is True:
                    exons[len(exons) - 1].query +=  sequence_match.groups()[0]
                else:
                    exons[len(exons) - 1].target += sequence_match.groups()[0]
                query_sequence_flag =               not query_sequence_flag
                continue
            
        return exons
        
    
if __name__ == '__main__':
    swParser = swAlignmentParser()
    swParser.setProteinFolder("ENSP00000311134")
    exons = swParser.parseOutput("/home/intern/Documents/sw_input/sw_output/ex6.swout", 1)
    for exon in sorted(exons, reverse=True):
        print exon.id, exon.match, exon.length, exon.score

        