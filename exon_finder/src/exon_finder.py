'''
Created on Mar 6, 2012

@author: marioot
'''
import os
import re
import sys
import ConfigParser
import csv
from os import path, system, remove
###############################
###############################
#  Tools configuration and setting paths
config_file         = "../../config.cfg"
config              = ConfigParser.RawConfigParser()
config.read(config_file)

blastn              = config.get('Blast cfg', 'blastn')
tblastn             = config.get('Blast cfg', 'tblastn')
swSharp             = config.get('SW#', 'swSharp')
sys.path.append(config.get('SW#', 'src'))

proteins_path       = ""
regions_path        = ""
blastout_path       = ""
swSharpout_path     = ""
exon_db             = ""
exon_individual_f   = ""
log_path            = ""
statistics_path     = ""
###############################
###############################
class Exon:
    score               = 0
    alignment_matches   = 0
    alignment_length    = 0
    alignment_start     = 0
    annotation          = ""
    def __init__(self, exon_id):
        self.id = id
###############################
###############################
#  Parses the description file
def parse_descr_file (descr_file) :
    species_name    ="";
    params          =[];
    all_species     = {};
    FILE            = open(descr_file,'r');
    while FILE :
        species_name = FILE.readline().rstrip('\n');
        if (not species_name):    # end of file
            break
        params_line = FILE.readline().rstrip('\n').split(' ')[1];
        params      = params_line.split(':');
        FILE.readline()
        all_species[species_name] = (species_name, params);
    FILE.close()
    return all_species
###############################
###############################
#  Species info to LOG file
def write_species_info (logfile, short_name, (species_name, params)):
    logfile.write("--SPECIES INFO:        {0} --\n\n".format(species_name))
    logfile.write("Source .fasta file:    {0}.fasta\n".format(regions_path + short_name))
    logfile.write("Type of sequence:      {0}\n".format(params[0]))
    logfile.write("Sequence ID:           {0}\n".format(params[2]))
    logfile.write("Region of interest:    {0}-{1}\n".format(params[3], params[4]))
    logfile.write("Strand:                {0}\n\n".format(params[5]))
###############################
###############################
#  Which exons have been found with blast tool
def parse_blastout(exons, blastout):
    results         = open(blastout, "r")
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
    results         = open(blastout, "r")
    for line in results.readlines():
        if(exon_counter == len(exons_found)):break;
        if(re.match('^>{0} '.format(exons_found[exon_counter]), line)):
            found_flag = True
        else:
            if(found_flag == True):
                reg = re.match(r'^\s*Score\s*=\s*(\d*.*|d*)\s*bits\s*\((\d*)\)', line)
                if(reg):
                    exons[exons_found[exon_counter]].score = int(reg.group(2))
                reg = re.match(r'^\s*Identities\s*=\s*(\d*)\/(\d*)', line)
                if(reg):
                    exons[exons_found[exon_counter]].alignment_matches  = int(reg.group(1))
                    exons[exons_found[exon_counter]].alignment_length   = int(reg.group(2))
                reg = re.match(r'^Query: (\d*)', line)
                if(reg):
                    exons[exons_found[exon_counter]].alignment_start    = int(reg.group(1))
                    found_flag = False
                    exon_counter += 1
    return [exons, exons_found]
###############################
###############################
#  Which exons have been found with swSharp tool
def parse_swSharpout(exons, swSharpout):
    #in development
    results     = open(swSharpout, "r")
    regex       = re.compile(r'Intervals: (\d+) (\d+).')
    
    for line in results:
        m = re.match(regex, line)
        if m is not None:
            lower = int(m.groups()[0])
            upper = int(m.groups()[1])
    return [lower, upper]
###############################
###############################
#  Align species DNA with exon database using blastn
def run_blastn (species_name) :
    blastout        = "{0}/dna/{1}.blastout".format(blastout_path, species_name)
    query_sequence  = "{0}/{1}.fa".format(regions_path, species_name)
    
    if (not path.isfile(query_sequence)) :
        print query_sequence, " file does not exist. "
        return -1;
    #blastn: species dna on exon database
    cmd = "{0} -d {1} -i {2} -o {3}".format(blastn, exon_db, query_sequence, blastout)
    system(cmd)
    return blastout
###############################
###############################
#  Align species protein with exon database using tblastn
def run_tblastn (species_name) :
    blastout        = "{0}/protein/{1}.blastout".format(blastout_path, species_name)
    query_sequence  = "{0}/{1}.fa".format(proteins_path, species_name)
    
    if (not path.isfile(query_sequence)) :
        print query_sequence, " file does not exist. "
        return -1;
    #tblastn: species dna on exon database
    cmd = "{0} -d {1} -i {2} -o {3}".format(tblastn, exon_db, query_sequence, blastout)
    system(cmd)
    return blastout
###############################
###############################
# Align species DNA with exons using CUDA Smith Waterman implementation
# Returns: {exon_num : [lower_bound, upper_bound], ...}
def run_SW (species_name, number_of_exons):
    # in development
    # todo: this function should return the output file of the tool, the database should contain all the exons
    swSharpout      = swSharpout_path + species_name + '.swSharpout'
    return swSharpout
    #####
    unimportant_file = "tmp.txt"
    query_sequence = regions_path + species_name + '.fa'
    sw_exons = {}
    for current_exon in range(1, number_of_exons + 1):
        exon_f = "{0}_{1}.fa".format(exon_individual_f, current_exon)
        swSharpout = swSharpout_path + "_{0}".format(current_exon)
        #pazi, promjena!!!
        cmd = "{0} -i {1} -j {2} --out {3} > {4}".format(swSharp, query_sequence, exon_f, swSharpout, unimportant_file)
        print cmd
        system(cmd)
        #sw_exons[current_exon] = parse_swSharpout_for_exons(swSharpout)
        remove(swSharpout)
    remove(unimportant_file)
    return sw_exons
###############################
###############################
# Given the exon information, infer about relative positions of exons and possible reasons why are some missing
# Returns a list of exons (given as a parameter) filled with appropriate annotation
def annotate_exons (exons_found, exons, sequence, number_of_exons):
    for current_exon in range(1, number_of_exons + 1):
        if current_exon in exons_found:
            exons[current_exon].annotation = "T"
        else:
            #find left edge
            left = 0
            reversed_range = range(0, current_exon+1)
            reversed_range.reverse() 
            for iterL in reversed_range:
                if iterL in exons_found:
                    left = exons[iterL].alignment_start + exons[iterL].alignment_length
                    break;
            #find right edge
            right = 0
            for iterR in range(current_exon, number_of_exons + 2):
                if iterR in exons_found:
                    right = exons[iterR].alignment_start
                    break;
            if(left != 0 and right != 0):
                matchObj = re.findall(r'N+', sequence[left:right])
                if len(matchObj) > 0 :
                    exons[current_exon].annotation = "N"
                else:
                    exons[current_exon].annotation = "?"
            else:
                if(left != 0):
                    matchObj = re.findall(r'N+', sequence[left:])
                    if len(matchObj) > 0 :
                        exons[current_exon].annotation = "ESN"
                    else:
                        exons[current_exon].annotation = "ES"
                else: 
                    if(right != 0):
                        matchObj = re.findall(r'N+', sequence[:right])
                        if len(matchObj) > 0 :
                            exons[current_exon].annotation = "ESN"
                        else:
                            exons[current_exon].annotation = "ES"
    return exons
###############################
###############################
#
def fetch_DNA_seq(species_name):
    fasta_f = regions_path + species_name + '.fa'
    DNA_f   = open (fasta_f, 'r');
    DNA_f.readline();
    seq     = DNA_f.read();
    seq     = "".join(seq.split('\n'))
    return seq
def initialize_exon_dictionary(number_of_exons):
    exons = {}
    for current_exon in range(1, number_of_exons + 1):
        exons[current_exon] = Exon(current_exon)
    return exons
###############################
###############################
#  
def find_exons (species_name, tool, number_of_exons):
    if tool == "tblastn":
        blastout                = run_tblastn(species_name)
        [exons, exons_found]    = parse_blastout(initialize_exon_dictionary(number_of_exons), blastout)
    elif tool == "blastn":
        blastout                = run_blastn(species_name)
        [exons, exons_found]    = parse_blastout(initialize_exon_dictionary(number_of_exons), blastout)
    else :
        swSharpout              = run_SW(species_name, number_of_exons)
        [exons, exons_found]    = parse_swSharpout(initialize_exon_dictionary(number_of_exons), swSharpout)
    sequence                    = fetch_DNA_seq(species_name)
    exons                       = annotate_exons(exons_found, exons, sequence, number_of_exons)
    return exons

def get_exon_lenghts(number_of_exons):
    exons       = {}
    exon_file   = open(exon_db, "r")
    for line in exon_file.readlines():
        m = re.match(r'>(\d+) exon length (\d+)', line)
        if m is not None:
            exons[int(m.groups()[0])] = int(m.groups()[1])
    return exons

def individual_statistic_annotation(exons, number_of_exons):
    # not in use ATM
    #[T, ?, N, ES, ESN]
    stat = {'T':0, '?':0, 'N':0, 'ES':0, 'ESN':0}
    print exons
    for current_exon in range(1, number_of_exons + 1):
        stat[exons[current_exon].annotation] = {'T': lambda: stat['T'] + 1,
                                                '?': lambda: stat['?'] + 1,
                                                'N': lambda: stat['N'] + 1,
                                                'ES': lambda: stat['ES'] + 1,
                                                'ESN': lambda: stat['ESN'] + 1}[exons[current_exon].annotation]()
    return stat

def format_individual_statistic_SW(exons_SW, number_of_exons):
    # not in use ATM
    stat = {}
    for current_exon in range(1, number_of_exons + 1):
        stat[current_exon] = 0.0
    return stat
    #
    stat = {}
    base_exon_sizes = get_exon_lenghts(number_of_exons)
    print "SW debug"
    print exons_SW
    print "base exon size:"
    print base_exon_sizes
    for current_exon in range(1, number_of_exons + 1):
        recognised_length = exons_SW[current_exon][1] - exons_SW[current_exon][0]
        print recognised_length
        stat[current_exon] = recognised_length / base_exon_sizes[current_exon]
    return stat
def generate_statistics_based_on_search(exons, base_exon_lengths):
    statistics = []
    for current_exon in range(1, len(exons) + 1):
        statistics.append([current_exon,
                           base_exon_lengths[current_exon],
                           exons[current_exon].score, 
                           exons[current_exon].alignment_matches, 
                           exons[current_exon].alignment_length,
                           exons[current_exon].annotation])
    return statistics
###############################
###############################
# Generates statistics 
# Returns: "<protein_id>, 
#           <species_name>, 
#           <type_of_search>, 
#           <exon>, 
#           <exon_length>, 
#           <score>, 
#           <number_of_matches>, 
#           <alignment_length>"
def generate_statistics(exons_via_proteins, exons_via_dna, exons_SW, all_species, protein_id):
    statistics_header = ["Protein_ID", 
                         "Species", 
                         "Type_of_search", 
                         "Exon_number", 
                         "Length", 
                         "Score", 
                         "Alignment_matches", 
                         "Alignment_length",
                         "Annotation"]
    tblastn_statistics = []
    blastn_statistics = []
    
    base_exon_length = get_exon_lenghts(len(exons_via_proteins))
    
    for species in sorted(all_species):
        tblastn_species_statistic = generate_statistics_based_on_search(exons_via_proteins[species], base_exon_length)
        blastn_species_statistic = generate_statistics_based_on_search(exons_via_dna[species], base_exon_length)
        for exon in tblastn_species_statistic:
            exon[1] /= 3    #for proteins
            exon.insert(0, "tblastn")
            exon.insert(0, species)
            exon.insert(0, protein_id)
        for exon in blastn_species_statistic:
            exon.insert(0, "blastn")
            exon.insert(0, species)
            exon.insert(0, protein_id)
        tblastn_statistics.append(tblastn_species_statistic)
        blastn_statistics.append(blastn_species_statistic)
        #print generate_statistics_based_on_search(exons_SW[species])
    ###WRITE TO CSV FILE###
    statout = csv.writer(open(statistics_path, 'wb+'), delimiter = ',')
    statout.writerow(statistics_header)
    for species in tblastn_statistics:
        for exon in species:
            statout.writerow(exon)
    for species in blastn_statistics:
        for exon in species:
            statout.writerow(exon)
    ####
    return [statistics_header, tblastn_statistics, blastn_statistics]

def generate_required_directories():
    if (not os.path.isdir(blastout_path)):
        os.makedirs("{0}/dna".format(blastout_path))
        os.makedirs("{0}/protein".format(blastout_path))

def main(argv=None):
    if argv is None:
        argv = sys.argv
    [number_of_exons, session_resource_path, protein_id] = [int(argv[1]), argv[2], argv[3]]
    #################
    #resulting proteins path
    global proteins_path
    global regions_path
    global blastout_path
    global swSharpout_path
    global exon_db
    global exon_individual_f
    global log_path
    global statistics_path
    proteins_path = "{0}/{1}".format(session_resource_path, config.get('Found proteins path', 'proteins'))
    regions_path = "{0}/{1}".format(session_resource_path, config.get('Gene regions path', 'regions'))
    blastout_path = "{0}/{1}".format(session_resource_path, config.get('Blastout path', 'blastout'))
    swSharpout_path = "{0}/{1}".format(session_resource_path, config.get('SW#', 'swSharpout'))
    exon_db = "{0}/{1}/exons.fa".format(session_resource_path, config.get('Exon database path', 'exons_path'))
    log_path = "{0}/{1}".format(session_resource_path, config.get('LOG', 'exon_finder'))
    statistics_path = "{0}/{1}.csv".format(session_resource_path, config.get('Statistics', 'exon_finder'))
    description_f = "{0}/{1}".format(session_resource_path, config.get('Session files', 'descr_output'))
    ################
    generate_required_directories()
    all_species = parse_descr_file(description_f)
    LOG = open(log_path, 'w')

    exons_via_proteins = {}
    exons_via_DNA={}
    exons_SW = {}
    for species in sorted(all_species):
        write_species_info(LOG, species, all_species[species])
        species_name= all_species[species][0];
        # 1)
        exons_via_proteins[species] = find_exons(species_name, "tblastn", number_of_exons)
        # 2)
        exons_via_DNA[species] = find_exons(species_name, "blastn", number_of_exons)
        # 3)
        #exons_SW[species] = find_exons(species_name, "SW", number_of_exons)
        
    print generate_statistics(exons_via_proteins, exons_via_DNA, exons_SW, all_species, protein_id)
    
    return 0
    
if __name__ == '__main__':
    sys.exit(main())