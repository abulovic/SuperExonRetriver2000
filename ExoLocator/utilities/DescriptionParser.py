'''
Created on Apr 18, 2012

@author: marioot
'''

# Python imports
import re

# utilities imports
from utilities.DirectoryCrawler import DirectoryCrawler


class DescriptionParser:
    '''
    Loads configuration files from the cfg directory
    '''
    def __init__(self):
        self.crawler = DirectoryCrawler()
    
    def get_gene_regions(self, protein_id):
        '''
        Parses the description file for the protein_id, and retrieves the information about protein locations for every species.
        Locations are stored as tuples (location_type, assembly, location_id, seq_begin, seq_end, strand)
        @param protein_id: protein_id for which protein ids of other species should be retrieved
        @return: (prot_ids_known, prot_ids_abinitio) - two dictionaries (key is species name, value is gene location data as described)
        '''
        (proteins_known, proteins_abinitio) = self.parse_descr_file(protein_id)
        
        genes_known = {}
        genes_abinitio = {}
        
        for key, value in proteins_known.items():
            genes_known[key] = list(value)[3:]
        for key, value in proteins_abinitio.items():
            genes_abinitio[key] = list(value)[1:]
            
        return genes_known, genes_abinitio

    def get_protein_ids(self, protein_id):
        '''
        Parses the description file for the protein_id, and retrieves only the protein ids for every species
        @param protein_id: protein_id for which protein ids of other species should be retrieved
        @return: (prot_ids_known, prot_ids_abinitio) - two dictionaries (key is species name, value is orthologous protein id)
        '''
        (proteins_known, proteins_abinitio) = self.parse_descr_file(protein_id)
            
        prot_ids_known = {}
        prot_ids_abinitio = {}
        
        for key, value in proteins_known.items():
            prot_ids_known[key] = list(value)[0]
        for key, value in proteins_abinitio.items():
            prot_ids_abinitio[key] = list(value)[0]
            
        return prot_ids_known, prot_ids_abinitio

    def get_species(self, protein_id):
        '''
        Parses the description file for the protein_id, and retrieves the list of species for which has reciprocal best search
        found a valid protein.
        @param protein_id: protein_id for which protein ids of other species should be retrieved
        @return: species_list - found by RBH
        '''
        (proteins_known, proteins_abinitio) = self.parse_descr_file(protein_id)
        
        species_list = proteins_known.keys()
        species_list.extend(proteins_abinitio.keys())
        
        for species in species_list:
            species = species.strip()
            
        return sorted(species_list)
    
    def get_separated_species (self, protein_id):
        (proteins_known, proteins_abinitio) = self.parse_descr_file(protein_id)
        species_known = proteins_known.keys()
        species_abinitio = proteins_abinitio.keys()
        return (species_known, species_abinitio)
    
    def get_strand_information (self, protein_id):
        '''
        Parses the description file for the protein_id and retrieves list of appropriate strands on which the proteins found by
        RBH are found.
        @param protein_id: protein_id for which protein ids of other species should be retrieved
        @return: strands: dictionary (species:strand) - found by RBH
        '''
        strands = {}
        (proteins_known, proteins_abinitio) = self.parse_descr_file(protein_id)
        for species, (spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand) in proteins_known.items():
            strands[species] = int(strand)
        for species, (spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand) in proteins_abinitio.items():
            strands[species] = int(strand)     
        return strands

    def parse_descr_file(self, protein_id):
        '''
        Function for parsing the description file associated with the protein_id.
        Description file contains two different types of entries: for known protein and for abinitio.
        Consequently, there are two formats that can be expected. They are both tab delimited.
        known_format:     species protein_id gene_id transcript_id location_type:assembly:location_id:seq_begin:seq_end:strand
        abinitio_format:  species protein_id location_type:assembly:location_id:seq_begin:seq_end:strand
        @param protein_id: protein for which the description file will be parsed. 
        @return: proteins_known_data - dictionary (key is species, value is a tuple of all the available data for that species protein)
                 abinitio_known_data (the same), two dictionaries are returned as a tuple
                 (spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
        @raise IOError: in case there is no description file present for the protein_id
        '''
        
        proteins_known_data = {}
        proteins_abinitio_data = {}
        
        descr_file_path = "{0}/{1}.descr".format(self.crawler.get_root_path(protein_id), protein_id)
        
        pattern_known = re.compile("(.*)\t(ENS.*)\t(.*)\t(.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
        pattern_abinitio = re.compile("(.*)\t(GEN.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
        
        try:
            descr_file = open(descr_file_path, 'r')
        except IOError:
            raise IOError("There is no description file present for protein: %s" % protein_id)
        
        for line in descr_file.readlines():
            line = line.strip()
            
            species_data = line.split()
            if len(species_data) == 5:
                species_name = species_data[0]
                spec_protein_id, gene_id, transcript_id = species_data[1:4]
                (location_type, assembly, location_id, seq_begin, seq_end, strand) = species_data[-1].split(":")
                proteins_known_data[species_name] = (spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)

                
            elif len(species_data) == 3:   
                species_name = species_data[0]
                (location_type, assembly, location_id, seq_begin, seq_end, strand) = species_data[-1].split(":")
                spec_protein_id = species_data[1]
                proteins_abinitio_data[species_name] = (spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
        
        descr_file.close()
        
        return proteins_known_data, proteins_abinitio_data
    
    def parse_description_file_general_info(self, protein_id):
        '''
        Similar as parse_description_file, but returns all the information. Required for generating a model structure.
        known: (species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
        abinitio: (species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)
        '''
        proteins_known_data = []
        proteins_abinitio_data = []
        
        descr_file_path = "{0}/{1}.descr".format(self.crawler.get_root_path(protein_id), protein_id)
        
        pattern_known = re.compile("(.*)\t(ENS.*)\t(.*)\t(.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
        pattern_abinitio = re.compile("(.*)\t(GEN.*)\t(.*):(.*):(.*):(.*):(.*):(.*)")
        
        try:
            descr_file = open(descr_file_path, 'r')
        except IOError:
            raise IOError("There is no description file present for protein: %s" % protein_id)
        
        for line in descr_file.readlines():
            line = line.strip()
            
            match = re.match(pattern_known, line)
            if match:
                (species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
                proteins_known_data.append((species_name, spec_protein_id, gene_id, transcript_id, location_type, assembly, location_id, seq_begin, seq_end, strand))
                
            match = re.match(pattern_abinitio, line)    
            if match:
                (species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand) = match.groups()
                proteins_abinitio_data.append((species_name, spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand))
        descr_file.close()
        
        return proteins_known_data, proteins_abinitio_data

def main():
    dp = DescriptionParser()
    print dp.get_species("ENSP00000298743")
    
if __name__ == '__main__':
    main()