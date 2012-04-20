'''
Created on Apr 15, 2012

@author: intern
'''
from subprocess                          import Popen, PIPE, STDOUT
from utilities.ConfigurationReader       import ConfigurationReader
from utilities.DescriptionParser         import DescriptionParser
from pipeline.utilities.DirectoryCrawler import DirectoryCrawler
import re, os, time
from utilities.Logger import Logger

def populate_sequence_exon_ensembl(protein_id):
    '''
    Populates the "/PROTEIN_ID/sequence/exon/ensembl/<species>.fa" 
    folder with fasta files containing a list of all the exons for
    a particular transcript. The data is aquired with a BioMartRemoteAccesss
    perl script.
    '''
    logger              = Logger.Instance()
    alignment_logger    = logger.get_logger('data_retrieval')
    perl_biomart_script = ConfigurationReader.Instance().get_value('data_retrieval', 'biomart_perl_script')
    exon_ensembl_path   = DirectoryCrawler().get_exon_ensembl_path(protein_id)
    tmp_xml_file        = "Query.xml"
    template_XML        = '''<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                    
        <Dataset name = "%s_gene_ensembl" interface = "default" >
            <Filter name = "ensembl_transcript_id" value = "%s"/>
            <Attribute name = "exon_chrom_start" />
            <Attribute name = "exon_chrom_end" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "gene_exon" />
            <Attribute name = "ensembl_exon_id" />
            <Attribute name = "strand" />
        </Dataset>
    </Query>
            '''
    try:
        (proteins_known, proteins_abinitio) = DescriptionParser().parse_descr_file(protein_id)
    except IOError, e:
        alignment_logger.error("{0}, {1}, , {2}".format(protein_id, 'ENSEMBL', e))
        return
    #(spec_protein_id, location_type, assembly, location_id, seq_begin, seq_end, strand)

    status_species_list = _read_failed_species(exon_ensembl_path)

    query_file_name     = "{0}/{1}".format(exon_ensembl_path, tmp_xml_file)
    failed_species_list = []
    for (species, data) in proteins_known.items():
        if species.strip() not in status_species_list and status_species_list != []:
            continue
        query_file      = open(query_file_name, 'w')
        transcript_id   = data[2]
        query = template_XML % (_query_name(species), transcript_id)
        query_file.write(query)
        query_file.close()
        exon_file_name  = "{0}/{1}.fa".format(exon_ensembl_path, species)
        biomart_command = "perl %s %s > %s" % (perl_biomart_script, query_file_name, exon_file_name)
        while (True):
            Popen(biomart_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            time.sleep(1)
            if (os.path.exists(exon_file_name)):
                break
        #LOGGING
        exon_file       = open(exon_file_name, 'r')
        exon_file_line  = exon_file.readline()
        exon_file.close()
        
        invalid_extension_pattern = re.compile("Query ERROR")
        if re.search(invalid_extension_pattern, exon_file_line) is None:
            continue
        alignment_logger.warning("{0}, {1}, {2}, {3}".format(protein_id, 'ENSEMBL', species.strip(), exon_file_line.strip()))
        os.remove(exon_file_name)
        failed_species_list.append(species.strip())
        
    try:
        os.remove(query_file_name)
    except OSError:
        pass
    
    if failed_species_list:
        _write_failed_species(exon_ensembl_path, failed_species_list)
        return False
    return True

def _write_failed_species(path, failed_species_list):
    '''
    Creates or updates the .status file for a respective sequence
    type with the names of species for which sequences haven't been found.
    '''
    status_file = open("{0}/.status".format(path), 'w')
    for species in failed_species_list:
        status_file.write("{0}\n".format(species))
    status_file.close()
    
    if failed_species_list == []:
        os.remove("{0}/.status".format(path))

def _read_failed_species(path):
    '''
    Reads the contents of the .status file and returns the list of species
    that have failed before.
    '''
    failed_species = []
    if os.path.exists("{0}/.status".format(path)):
        status_file = open("{0}/.status".format(path), 'r')
        for line in status_file.readlines():
            failed_species.append(line.strip())
        status_file.close()
    return failed_species  
  
def _query_name(species_name):
    '''
    Helper function that derives a name of the species for the XML query.
    '''
    name_parts = re.split('_', species_name)
    return "{0}{1}".format(name_parts[0][0].lower(), name_parts[1])

def main ():
    populate_sequence_exon_ensembl("ENSP00000311134")

if __name__ == '__main__':
    main()