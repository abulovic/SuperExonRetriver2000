'''
Created on Apr 23, 2012

@author: intern
'''
import re
from utilities.ConfigurationReader import ConfigurationReader

def extract_failed_protein_description(output_file_name = None):
    cr = ConfigurationReader.Instance()
    
    failed_protein_file_path = cr.get_value('input', 'failed_proteins')
    protein_descr_file_path = cr.get_value('input', 'protein_description')
    
    try:
        failed_protein_file = open(failed_protein_file_path, 'r')
    except IOError:
        raise IOError ("Value of failed_proteins, section input in the directory_tree.cfg is pointing to a non-existant file.")
    try:
        protein_descr_file = open(protein_descr_file_path, 'r')
    except IOError:
        raise IOError ("Value of protein_description, section input in the directory_tree.cfg is pointing to a non-existant file.")
    
    # read the protein ids from failed protein file
    failed_protein_list = list (prot_id.strip() for prot_id in failed_protein_file.readlines())
    
    print "Unable to retrieve %d proteins" % len(failed_protein_list)
    
    pattern = re.compile("\d+\s+(ENSP\w*\d+)\s+(\S.*)")
    protein_descr_dict = {}
    uncharacterized = {}
    immune_system = {}
    olfactory = {}
    others = {}
    pseudo = {}
    no_descr = {}
    for line in protein_descr_file.readlines():
        line = line.strip()
        #print line
        protein_descr_match = re.match(pattern, line)
        if protein_descr_match:
            (prot_id, descr) = protein_descr_match.groups()
            protein_descr_dict[prot_id] = descr          
    
    for failed_prot_id in failed_protein_list:
        descr = protein_descr_dict.get(failed_prot_id, "No description available")
        
        if "uncharacterized" in descr.lower() or "open reading frame" in descr.lower():
            uncharacterized[failed_prot_id] = descr
            continue
        if "olfactory" in descr.lower():
            olfactory[failed_prot_id] = descr
            continue
        if "t cell" in descr.lower() or "immunoglobulin" in descr.lower() or "histocompatibility complex" in descr.lower():
            immune_system[failed_prot_id] = descr
            continue
        if "pseudogene" in descr.lower():
            pseudo[failed_prot_id] = descr
            continue
        if "No description available" in descr:
            no_descr[failed_prot_id] = descr
            continue
        others[failed_prot_id] = descr
        
    for (pid,descr) in others.items():
        print pid, descr 
        
    protein_descr_file.close()
    failed_protein_file.close()
        
    
    
    
if __name__ == '__main__':
    extract_failed_protein_description()
