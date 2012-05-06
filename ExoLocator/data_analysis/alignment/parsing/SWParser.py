'''
Created on May 1, 2012

@author: intern
'''

import re, os

from utilities.Logger                       import Logger
from data_analysis.base.Exon                import Exon
from pipeline.utilities.DirectoryCrawler    import DirectoryCrawler


def parse_SW_output (ref_protein_id, species, sw_type):
    
    logger              = Logger.Instance()
    containers_logger   = logger.get_logger('containers')
    dc                  = DirectoryCrawler()
    
    if sw_type.lower() == "sw_gene":
        swout_file_path = dc.get_SW_gene_path(ref_protein_id)
    elif sw_type.lower() == "sw_exon":
        swout_file_path = dc.get_SW_exon_path(ref_protein_id)
    else:
        raise KeyError ("There is no known swout path for type %s" % sw_type)
    swout_file_path += "/%s.swout" % species
    
    if not os.path.isfile(swout_file_path):
        containers_logger.error ("{0}, {1}, {2}, no swout file".format(ref_protein_id, species, sw_type))
        return False
    
    swout_file = open(swout_file_path, 'r')
    
    # status boolean variables
    parsing_query_seq = True
    
    # patterns for matching
    header_pattern      = re.compile ("Name: >(\d+)\|(\d+)\|(ENS\w+)\|(ENS\w+)\|([-]*1)")
    #Intervals: 1207047 1207087 30 69 (+) strand 
    intervals_pattern   = re.compile ("Intervals:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\([+-]\)\s+strand")
    #Identity: 31/41 (75.6%)
    identity_pattern    = re.compile ("Identity:\s+(\d+)/(\d+).*")
    #Similarity: 40/41 (97.6%)
    similarity_pattern  = re.compile ("Similarity:\s+(\d+)/.*")
    #Gaps: 1/41 (2.4%)
    gaps_pattern        = re.compile ("Gaps:\s+(\d+)/\d+.*")
    # sequence pattern
    sequence_pattern    = re.compile ("\s*(\d+)\s+([ATCG-]+)\s+(\d+).*")
    
    exon_dict = {}
    ref_exon_id     = ""
    identities      = 0
    positives       = 0
    gaps            = 0
    sbjct_start     = 0
    sbjct_end       = 0
    query_start     = 0
    query_end       = 0
    length          = 0
    query_sequence  = ""
    sbjct_sequence  = ""
    exon = Exon(sw_type, "")
    
    for line in swout_file.readlines():
        line = line.strip()
        header_match = re.match(header_pattern, line)
        if header_match:
            #create new exon and start a new one
            if ref_exon_id:
                exon.set_alignment_info(int(identities), 
                                        int(positives), 
                                        int(gaps), 
                                        int(sbjct_start), 
                                        int(sbjct_end), 
                                        int(query_start), 
                                        int(query_end), 
                                        int(length), 
                                        sbjct_sequence,
                                        query_sequence)
                if ref_exon_id in exon_dict:
                    exon_dict[ref_exon_id].append(exon)
                else:
                    exon_dict[ref_exon_id] = [exon]
                 
            
            ref_exon_id = header_match.groups()[3]
            exon = Exon(sw_type, ref_exon_id)
            sequence = ""
            
        intervals_match = re.match (intervals_pattern, line)
        if intervals_match:
            (query_start, query_end, sbjct_start, sbjct_end) = intervals_match.groups()
            
        identity_match = re.match (identity_pattern, line)
        if identity_match:
            (identities, length) = identity_match.groups()
            
        similarity_match = re.match(similarity_pattern, line)
        if similarity_match:
            positives = similarity_match.groups()[0]
            
        gaps_match = re.match (gaps_pattern, line)
        if gaps_match:
            gaps = gaps_match.groups()[0]
            
        sequence_match = re.match (sequence_pattern, line)
        if sequence_match:
            if parsing_query_seq:
                query_sequence += sequence_match.groups()[1].strip()
                parsing_query_seq = False
            else:
                sbjct_sequence += sequence_match.groups()[1].strip()
                parsing_query_seq = True
                
    exon.set_alignment_info(int(identities), 
                            int(positives), 
                            int(gaps), 
                            int(sbjct_start), 
                            int(sbjct_end), 
                            int(query_start), 
                            int(query_end), 
                            int(length), 
                            sbjct_sequence,
                            query_sequence)
    if query_sequence:
        if ref_exon_id in exon_dict:
            exon_dict[ref_exon_id].append(exon)
        else:
            exon_dict[ref_exon_id] = [exon]
                
                
    return exon_dict
    
    
    
    
if __name__ == '__main__':
    exon_dict = parse_SW_output ("ENSP00000365108", "Homo_sapiens", "sw_exon")
    for ref_exon_id, exon_list in exon_dict.items():
        print ref_exon_id
        for al_exon in exon_list:
            print "\tSbjct_start:  %d" % al_exon.alignment_info["sbjct_start"]
            print "\tQuery_start:  %d" % al_exon.alignment_info["query_start"]
            print "\tAlign_length: %d" % al_exon.alignment_info["length"]
            print