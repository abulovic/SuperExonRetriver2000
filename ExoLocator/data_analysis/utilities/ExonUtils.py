'''
Created on Jun 7, 2012

@author: intern
'''
from data_analysis.containers.ProteinContainer import ProteinContainer
from data_analysis.containers.DataMapContainer import DataMapContainer
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from data_analysis.base.EnsemblExon import EnsemblExon

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
    except KeyError:
        return None
    
    new_exons =[]
    
    for exon in exons:
       
        new_exon = EnsemblExon((exon.ref_protein_id, exon.species), exon.exon_id, exon.start, exon.stop, exon.strand, exon.sequence)
        new_exon.set_exon_ordinal(exon.ordinal)
        
    
    return new_exons
        