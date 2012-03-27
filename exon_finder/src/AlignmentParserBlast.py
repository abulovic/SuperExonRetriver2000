'''
Created on Mar 27, 2012

@author: Mario, anana
'''

import re
from AlignmentParser import AlignmentParser

class blastAlignmentParser(AlignmentParser):
    '''
    Class offers methods for parsing various blast outputs.
    Refers to blastn / tblastn.
    '''
    
    def __init__(self, algorithmName="blast"):
        '''
        @param algorithmName: optional - algorithm name, blast by default
        '''
        self.algorithmName = algorithmName
        super(blastAlignmentParser, self).__init__()
    
        
        
    def parseOutput (self, alignmentOutputFile, numberOfExons):
        '''
        Parses blastn / tblastn output
        @param alignmentOutputFile: the result of blastn (tblastn) alignment
        @param numberOfExons: number of exons in reference species
        @return: [exons, exons_found] First is the list of all the exons
                (Exon objects), the second a list of number IDs of found exons.
        '''
        
        exons = self.initializeExonDictionary(numberOfExons)
        results         = open(alignmentOutputFile, "r")
        reading_flag    = False
        # Identify exons found by blast
        out             = "";
        exons_found     = [];
        
        for line in results.readlines():
            if(re.match('^Sequences producing significant alignments:', line) or re.match('^>', line)):
                if(reading_flag) :
                    break
                reading_flag = True
                out = ""
            else:
                out += line
                
        for line in out.split("\n"):
            if (line != ""):
                exon_id = int(line.split(" ")[0])
                exons_found.append(exon_id)
                
        # Get information about those exons
        exon_counter    = 0
        found_flag      = False
        results         = open(alignmentOutputFile, "r")
        for line in results.readlines():
            if(exon_counter == len(exons_found)):break;
            if(re.match('^>{0} '.format(exons_found[exon_counter]), line)):
                found_flag = True
            else:
                if(found_flag == True):
                    reg = re.match(r'^\s*Score\s*=\s*(\d*.*|d*)\s*bits\s*\((\d*)\)', line)
                    if(reg):
                        exons[exons_found[exon_counter]].set_score(int(reg.group(2)))
                        #exons[exons_found[exon_counter]].score = int(reg.group(2))
                    reg = re.match(r'^\s*Identities\s*=\s*(\d*)\/(\d*)', line)
                    if(reg):
                        exons[exons_found[exon_counter]].set_identity (int(reg.group(1)), int(reg.group(2)))
                        #exons[exons_found[exon_counter]].alignment_matches  = int(reg.group(1))
                        #exons[exons_found[exon_counter]].alignment_length   = int(reg.group(2))
                    reg = re.match(r'^Query: (\d*)', line)
                    if(reg):
                        exons[exons_found[exon_counter]].alignment_start    = int(reg.group(1))
                        found_flag = False
                        exon_counter += 1
        return [exons, exons_found]
        
        
if __name__ == '__main__':
    blastParser = blastAlignmentParser()
    blastParser.setProteinFolder("ENSP00000311134")
    [exons, exons_found] = blastParser.parseOutput("/home/intern/Project/workspaceBIO/sessions/ENSP00000311134/blastout/dna/Ailuropoda_melanoleuca.blastout", 15)
    for key in exons_found:
        exon = exons[key]
        print exon
        print exon.match, exon.length
        
        
        
        
        