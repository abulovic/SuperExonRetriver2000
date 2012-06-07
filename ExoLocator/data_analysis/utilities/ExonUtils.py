'''
Created on Jun 7, 2012

@author: intern
'''
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer

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
        exon_seq = exon.sequence
        
        exon_translated = False
        exon_start = False
        exon_stop = False
        
        for frame in (0,1,2):
            translation = str(exon_seq[frame:].translate())
            longest_translation = LongestCommonSubstring(protein_sequence, translation)
            # exon u sredini
            if longest_translation == translation:
                exon_translated = True
                break
                
            if protein_sequence.startswith(longest_translation) and translation.endswith(longest_translation):
                exon_start = True
                break
       
            if protein_sequence.endswith(longest_translation) and translation.startswith(longest_translation):
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
        