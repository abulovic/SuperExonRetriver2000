'''
Created on Apr 20, 2012

@author: marioot
'''
from utilities                  import FileUtilities
from pipeline.alignment import Alignments
import sys

def populate_referenced_species_databases(protein_list, referenced_species):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['ENSEMBL_EXON_RETRIEVAL'] == 'FAILED':
            print "ABORTING {0} DATABASE FORMATTING: ENSEMBL_EXON_RETRIEVAL has a FAILED status!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'REF_SP_DB_FORMATTING', 'FAILED')
            continue
        try:
            if FileUtilities.read_status_file(protein_id)['REF_SP_DB_FORMATTING'] == 'OK':
                print "SKIPPING {0} DATABASE FORMATTING: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "DATABASE FORMATTING: {0} for reference species {1}".format(protein_id, referenced_species)
        if Alignments.generate_referenced_species_database(protein_id, referenced_species):
            FileUtilities.update_entry_in_status_file(protein_id, 'REF_SP_DB_FORMATTING', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'REF_SP_DB_FORMATTING', 'FAILED') 

def populate_blastn_alignments(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['EXP_GENE_RETRIEVAL'] == 'FAILED' or FileUtilities.read_status_file(protein_id)['REF_SP_DB_FORMATTING'] == 'FAILED':
            print "ABORTING {0} BLASTN: some resources have FAILED stats!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'BLASTN_ALIGNMENT', 'FAILED')
            continue
        try:
            if FileUtilities.read_status_file(protein_id)['BLASTN_ALIGNMENT'] == 'OK':
                print "SKIPPING {0} BLASTN: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "ALIGNING BLASTN: {0}".format(protein_id)
        if Alignments.generate_blastn_alignments(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'BLASTN_ALIGNMENT', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'BLASTN_ALIGNMENT', 'PARTIAL') 
      
def populate_tblastn_alignments(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['PROTEIN_RETRIEVAL'] == 'FAILED' or FileUtilities.read_status_file(protein_id)['REF_SP_DB_FORMATTING'] == 'FAILED':
            print "ABORTING {0} TBLASTN: some resources have FAILED stats!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'TBLASTN_ALIGNMENT', 'FAILED')
            continue
        try:
            if FileUtilities.read_status_file(protein_id)['TBLASTN_ALIGNMENT'] == 'OK':
                print "SKIPPING {0} TBLASTN: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "ALIGNING TBLASTN: {0}".format(protein_id)
        if Alignments.generate_tblastn_alignments(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'TBLASTN_ALIGNMENT', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'TBLASTN_ALIGNMENT', 'PARTIAL') 

def populate_SW_gene_alignments(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['EXP_GENE_RETRIEVAL'] == 'FAILED' or FileUtilities.read_status_file(protein_id)['REF_SP_DB_FORMATTING'] == 'FAILED':
            print "ABORTING {0} SW_GENE: some resources have FAILED stats!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_GENE_ALIGNMENT', 'FAILED')
            continue
        try:
            if FileUtilities.read_status_file(protein_id)['SW_GENE_ALIGNMENT'] == 'OK':
                print "SKIPPING {0} SW_GENE: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "ALIGNING SW_GENE: {0}".format(protein_id)
        if Alignments.generate_SW_gene_alignments(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_GENE_ALIGNMENT', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_GENE_ALIGNMENT', 'PARTIAL') 

def populate_SW_exon_alignments(protein_list):
    for protein_id in protein_list:
        if FileUtilities.read_status_file(protein_id)['ENSEMBL_EXON_RETRIEVAL'] == 'FAILED' or FileUtilities.read_status_file(protein_id)['REF_SP_DB_FORMATTING'] == 'FAILED':
            print "ABORTING {0} SW_EXON: some resources have FAILED stats!".format(protein_id)
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_EXON_ALIGNMENT', 'FAILED')
            continue
        try:
            if FileUtilities.read_status_file(protein_id)['SW_EXON_ALIGNMENT'] == 'OK':
                print "SKIPPING {0} SW_EXONs: .status file -> OK!".format(protein_id)
                continue
        except KeyError:
            pass
        print "ALIGNING SW_EXON: {0}".format(protein_id)
        if Alignments.generate_SW_exon_alignments2(protein_id):
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_EXON_ALIGNMENT', 'OK')
        else:
            FileUtilities.update_entry_in_status_file(protein_id, 'SW_EXON_ALIGNMENT', 'PARTIAL') 

def main ():
    referenced_species = "Homo_sapiens"
    
    protein_list_raw = FileUtilities.get_protein_list()
    protein_list = []
    for protein_tuple in protein_list_raw:
        protein_list.append(protein_tuple[0])
    
    if(len(sys.argv) < 1):
        print "Usage: {0} <blastn | tblastn | SW_gene | SW_exon | all> \n".format(sys.argv[0])
        exit
    mode = sys.argv[1]
    
    populate_referenced_species_databases(protein_list, referenced_species)
    
    if (mode == "blastn"):
        populate_blastn_alignments(protein_list)
    elif (mode == "tblastn"):
        populate_tblastn_alignments(protein_list)
    elif (mode == "SW_gene"):
        populate_SW_gene_alignments(protein_list)
    elif (mode == "SW_exon"):
        populate_SW_exon_alignments(protein_list)
    elif (mode == "all"):
        populate_blastn_alignments(protein_list)
        populate_tblastn_alignments(protein_list)
        populate_SW_gene_alignments(protein_list)
        populate_SW_exon_alignments(protein_list)
    else:
        print "Usage: {0} <blastn | tblastn | SW_gene | SW_exon | all> \n".format(sys.argv[0])
        exit
    
if __name__ == '__main__':
    main()