'''
Created on May 1, 2012

@author: intern
'''
# Python imports
import re, os

# utilities imports
from utilities.Logger                       import Logger
from utilities.DirectoryCrawler             import DirectoryCrawler

# data analysis imports
from data_analysis.base.Exon                import Exon


def parse_SW_output (ref_protein_id, species, sw_type):
    '''
    Parses the output from the SW# command line application.
    (suitable for version as it was distributed on May 1st, 2012)
    
    @param sw_type: sw_exon/sw_gene
    @return: dictionary of alignment exons. The keys are referent exon IDs, and 
    values are lists of all the alignment exons which correspond to the certain
    reference exon 
    '''
    
    logger              = Logger.Instance()
    containers_logger   = logger.get_logger('containers')
    dc                  = DirectoryCrawler()
    
    # determine the swout file path
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
    sequence_pattern    = re.compile ("\s*(\d+)\s+([ATCGN-]+)\s+(\d+).*")
    
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
    exon = Exon(sw_type, "", ref_protein_id, species)
    
    for line in swout_file.readlines():
        
        line = line.strip()
        header_match = re.match(header_pattern, line)
        if header_match:
            #add the current exon and start a new one
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
            exon = Exon(sw_type, ref_exon_id, ref_protein_id, species)
            parsing_query_seq = True
            query_sequence = ""
            sbjct_sequence = ""
        
        # intervals    
        intervals_match = re.match (intervals_pattern, line)
        if intervals_match:
            (query_start, query_end, sbjct_start, sbjct_end) = intervals_match.groups()
            
        # identities
        identity_match = re.match (identity_pattern, line)
        if identity_match:
            (identities, length) = identity_match.groups()
            
        # similarities
        similarity_match = re.match(similarity_pattern, line)
        if similarity_match:
            positives = similarity_match.groups()[0]
            
        # gaps
        gaps_match = re.match (gaps_pattern, line)
        if gaps_match:
            gaps = gaps_match.groups()[0]
            
        # sequence
        sequence_match = re.match (sequence_pattern, line)
        if sequence_match:
            sequence_to_append = sequence_match.groups()[1].strip()
            if parsing_query_seq:
                query_sequence += sequence_to_append
                parsing_query_seq = False
            else:
                sbjct_sequence += sequence_to_append
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
    exon_dict = parse_SW_output ("ENSP00000253108", "Ailuropoda_melanoleuca", "sw_gene")
    for ref_exon_id, exon_list in exon_dict.items():
        print ref_exon_id
        for al_exon in exon_list:
            print "\tSbjct_start:  %d" % al_exon.alignment_info["sbjct_start"]
            print "\tQuery_start:  %d" % al_exon.alignment_info["query_start"]
            print "\tAlign_length: %d" % al_exon.alignment_info["length"]
            print "\tSbjct seq:    %s" % al_exon.alignment_info["sbjct_seq"]
            print "\tQuery seq:    %s" % al_exon.alignment_info["query_seq"]
            
            print