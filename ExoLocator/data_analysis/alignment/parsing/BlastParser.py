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
        return False
        
    file_handle = open(blast_file, 'r')
    
    # parse blastn output  
    blastn_record = NCBIXML.read(file_handle)
    
    exon_dict = {}
    exon_pattern = re.compile("(\d+)\|(\d+)\|(ENS\w+)\|(ENS\w+)\|([+-]1)")
    
    for alignment in blastn_record.alignments:
        (blast_info, exon_info) = alignment.title.split()
        pattern_match = re.match(exon_pattern, exon_info)
        ref_exon_id = pattern_match.groups()[3]
        for hsp in alignment.hsps:
            exon = Exon(blast, ref_exon_id)
            if type(hsp.gaps) is int:
                gaps = hsp.gaps
            elif type(hsp.gaps) is tuple:
                if not hsp.gaps[0]:
                    gaps = 0
            exon.set_alignment_info (hsp.identities, 
                                                  hsp.positives, 
                                                  gaps, 
                                                  hsp.sbjct_start, 
                                                  hsp.sbjct_start + len(hsp.sbjct) -1,
                                                  hsp.query_start,
                                                  hsp.query_start + len(hsp.sbjct) -1,
                                                  len(hsp.sbjct),
                                                  hsp.sbjct)
            if not ref_exon_id in exon_dict:
                exon_dict[ref_exon_id] = [exon]
            else:
                exon_dict[ref_exon_id].append(exon)
        
    file_handle.close()
    return exon_dict



def main():
    exon_dict = parse_blast_output ("ENSP00000341765", "Ailuropoda_melanoleuca", "blastn")
    for (ref_exon_id, exon) in exon_dict.items():
        print ref_exon_id
        print "NUM OF ALIGNMENTS: %d" % len(exon.reference_exons[ref_exon_id])
        for alignment_info in exon.reference_exons[ref_exon_id]:
            print "\tAlignment coord: %d-%d" % (alignment_info["start"], alignment_info["stop"])
            print "\tAlignment ids:   %d/%d" % (alignment_info["identities"], alignment_info["ref_exon_len"])
            print "\tAlignment seq:   %s" % alignment_info["sequence"]
            print
        

if __name__ == '__main__':
    main()
