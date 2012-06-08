'''
Created on Jun 7, 2012

@author: intern
'''
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def LongestCommonSubstring(S1, S2):
    M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))]
    longest, x_longest = 0, 0
    for x in xrange(1,1+len(S1)):
        for y in xrange(1,1+len(S2)):
            if S1[x-1] == S2[y-1]:
                M[x][y] = M[x-1][y-1] + 1
                if M[x][y]>longest:
                    longest = M[x][y]
                    x_longest  = x
            else:
                M[x][y] = 0
    return S1[x_longest-longest: x_longest]

def longest_border_substring (s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    
    if l1 >= l2:
        longer_str = s1
        shorter_str = s2
    else:
        longer_str = s2
        shorter_str = s1
    
    longest_substring = ""
    
    # provjeri pocetak
    start_index = 0
    end_index = len(shorter_str)
    
    while (True):
        if start_index == len(shorter_str)-1:
            break
        substring = shorter_str [start_index:]
        if longer_str.startswith(substring):
            longest_substring = substring
            break
        start_index += 1
    
    # provjeri kraj
    
    while (True):
        if end_index == 0:
            break
        substring = shorter_str[0:end_index]
        if longer_str.endswith(substring):
            if len(substring) > len(longest_substring):
                longest_substring = substring
            break
        end_index -= 1
        
    return longest_substring
        
    

def remove_UTR_ensembl_exons (protein_id, species, exons):
    
    '''
    Removes the untranslated regions from ensembl exons
    This is for the purpose of statistics generating
    '''
    pc  = ProteinContainer.Instance()
    dmc = DataMapContainer.Instance()
    dm_key = (protein_id, species)
    try:
        dm = dmc.get(dm_key)
        prot = pc.get(dm.protein_id)
    except KeyError:
        return None
    
    new_exons =[]
    protein_sequence = str(prot.get_sequence_record().seq)
    
    for exon in exons:
        exon_seq = Seq(str(exon.sequence), IUPAC.ambiguous_dna)

        exon_translated = False
        exon_start = False
        exon_stop = False
        
        for frame in (0,1,2):
            translation = str(exon_seq[frame:].translate())
            longest_translation = LongestCommonSubstring(protein_sequence, translation)
            edge_translation = longest_border_substring(protein_sequence, translation)
            # exon u sredini
            if longest_translation == translation:
                exon_translated = True
                break
            
            elif protein_sequence.startswith(longest_translation) and translation.endswith(longest_translation):
                exon_start = True
                break
       
            elif protein_sequence.endswith(longest_translation) and translation.startswith(longest_translation):
                exon_stop = True
                break
    
            elif protein_sequence.startswith(edge_translation) and translation.endswith(edge_translation):
                exon_start = True
                break
       
            elif protein_sequence.endswith(edge_translation) and translation.startswith(edge_translation):
                exon_stop = True
                break
            
            
            
        new_exon = exon
        
        if exon_start:
            new_exon.sequence = exon_seq[(len(translation) - len(longest_translation))*3 + frame:]
            new_exon.start = exon.start + frame + len(longest_translation)*3
        elif exon_translated:
            new_exon.sequence = exon_seq
        elif exon_stop:
            new_exon.sequence = exon_seq[:frame+len(longest_translation)*3]
            new_exon.stop = exon.start + frame + len(longest_translation)*3
        else:
            continue
        
        new_exons.append(new_exon)
    
    return new_exons
        