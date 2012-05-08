'''
Created on May 1, 2012

@author: intern
'''
import re, os

from Bio.Blast                              import NCBIXML

from utilities.Logger                       import Logger
from data_analysis.base.Exon                import Exon
from pipeline.utilities.DirectoryCrawler    import DirectoryCrawler


def parse_blast_output (ref_protein_id, species, blast):
    '''
    Returns a dictionary. Key is reference species exon_id, and 
    the value is list of corresponding alignments
    '''
    
    logger              = Logger.Instance()
    containers_logger   = logger.get_logger('containers')
    dc                  = DirectoryCrawler()
    
    if blast == "blastn":
        blast_file = "{0}/{1}.blastout".format(dc.get_blastn_path(ref_protein_id), species)
    else:
        blast_file = "{0}/{1}.blastout".format(dc.get_tblastn_path(ref_protein_id), species)
        
    if not os.path.isfile(blast_file):
        containers_logger.error ("{0}, {1}, {2}, no blastout file".format(ref_protein_id, species, blast))
        return None
        
    file_handle = open(blast_file, 'r')
    
    # parse blastn output  
    try:
        blastn_record = NCBIXML.read(file_handle)
    except ValueError:
        containers_logger.error("%s,%s,%s,No hits found" % (ref_protein_id, species, blast))
        return None
    
    exon_dict = {}
    exon_pattern = re.compile("(\d+)\|(\d+)\|(ENS\w+)\|(ENS\w+)\|([-]*1)")
    
    for alignment in blastn_record.alignments:
        (blast_info, exon_info) = alignment.title.split()
        pattern_match = re.match(exon_pattern, exon_info)
        ref_exon_id = pattern_match.groups()[3]
        
        # limit alignments to 10 hsps
        
        num_of_hsps = 0
        
        for hsp in alignment.hsps:
            # limit!
            if blast == "blastn":
                (query_frame, hit_frame) = hsp.frame
                if query_frame == -1 or hit_frame == -1:
                    continue 
            if num_of_hsps == 5:
                break
            num_of_hsps += 1
            
            exon = Exon(blast, ref_exon_id)
            if type(hsp.gaps) is int:
                gaps = hsp.gaps
            elif type(hsp.gaps) is tuple:
                if not hsp.gaps[0]:
                    gaps = 0
            exon.set_alignment_info ( hsp.identities, 
                                      hsp.positives, 
                                      gaps, 
                                      hsp.sbjct_start, 
                                      hsp.sbjct_start + len(hsp.sbjct) -1,
                                      hsp.query_start,
                                      hsp.query_start + len(hsp.sbjct) -1,
                                      len(hsp.sbjct),
                                      hsp.sbjct,
                                      hsp.query)
            if not ref_exon_id in exon_dict:
                exon_dict[ref_exon_id] = [exon]
            else:
                exon_dict[ref_exon_id].append(exon)
        
    file_handle.close()
    return exon_dict



def main():
    exon_dict = parse_blast_output ("ENSP00000331363", "Callithrix_jacchus", "blastn")
    for (ref_exon_id, al_exons) in exon_dict.items():
        print ref_exon_id
        for al_exon in al_exons:
            print al_exon.alignment_info["identities"]
            print al_exon.alignment_info["length"]
        
        

if __name__ == '__main__':
    main()
