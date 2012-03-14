'''
Created on Dec 8, 2011

@author: Mario
'''

from os import system
import re
from os import path

# path to descr file
descr_file = "/Users/Mario/Desktop/project/ExonFinder/ensembl_mutual_best_parts/src/mammals.descr";

# path to perl file that can grab single animal slice i do direktorija s regijama
grab_slice = "../../ensembl_search/grab_single_slice.pl";
regions_path = "/Users/Mario/Desktop/project/ExonFinder/ensembl_search/Regions/"
proteins_path = "./proteins/"

# number of exons
number_of_exons = 15;

# path to log file
log_file = "./exon_finder.log"


blastn = "/Users/Mario/Desktop/Bioinf/blastn-2.2.16/bin/blastall -p blastn -e 1.e-2"
tblastn = "/Users/Mario/Desktop/Bioinf/blastn-2.2.16/bin/blastall -p tblastn -e 1.e-2"
exon_db = "/Users/Mario/Desktop/project/ExonFinder/ExonFinder/src/human_parf_exons.fasta"

#ovako bi trebao izgledati input file
query_sequence = "MYO_LUC_scaffold.fasta"
outfile = "tmp.blastout"

#statistika
#[T, N, X, NX, U, valid_species]
statistics = [0 for i in range(6)]
statistics_w = [0 for i in range(6)]
exon_statistics = [0 for i in range(number_of_exons)]
new_exon_discovery_statistics = [0 for i in range(number_of_exons)]
protein_exon_discovery_statistics = [0 for i in range(number_of_exons)]

new_exon_discovery_species = {}
protein_exon_discovery_species = {}
#Baza se prije koristenja treba formatirati i 
#smjestiti sve elemente baze u isti folder gdje se nalazi i query sequence
#Naredba za formatiranje baze: <staza do formatdb> -i <fasta datoteka> -p F -o T

def test_for_N_regions (fasta_file) :
	if (not path.isfile(fasta_file)) :
		return -1
		
	file = open (fasta_file, 'r');
	file.readline();
	seq = file.read();
	seq = "".join(seq.split('\n'))
	matchObj = re.findall(r'N+', seq)
	if len(matchObj) > 0 :
		print "There are ", len(matchObj), " sequences of Ns present. Some data might be missing."
	return len(matchObj)

def find_exons (params, species_name, short_name) :
	outfile = './blastout/'+short_name+'.blastout'
	query_sequence = regions_path + short_name + '.fasta';
	
	if (not path.isfile(query_sequence)) :
		print query_sequence, " file does not exist. "
		return -1;
	
	
	print "Current query being executed on ", query_sequence, "."
	
	#blastn exons against found scaffolds in other species
	cmd = "%s -d %s -i %s -o %s" % (blastn, exon_db, query_sequence, outfile)
	print cmd
	system(cmd)
	
	results = open(outfile, "r")
	reading = False
	
	out = "";
	exons = [];
	for line in results.readlines():
		if(re.match('^Sequences producing significant alignments:', line) or re.match('^>', line)):
			if(reading) :
				break
			reading = True
			out = "";
		else:
			out += line

	for line in out.split("\n"):
		if (line != ""):
			exons.append(line.split(" ")[0])
	return exons
	
def find_exons_prot (params, species_name, short_name) :
	outfile = './blastout/proteins/'+short_name+'.blastout'
	query_sequence = proteins_path + short_name + '.fasta';
	
	if (not path.isfile(query_sequence)) :
		print query_sequence, " file does not exist. "
		return -1;
	
	
	print "Current PROTEIN query being executed on ", query_sequence, "."
	
	#blastn exons against found proteins in other species
	cmd = "%s -d %s -i %s -o %s" % (tblastn, exon_db, query_sequence, outfile)
	print cmd
	system(cmd)
	
	results = open(outfile, "r")
	reading = False
	
	out = "";
	exons = [];
	for line in results.readlines():
		if(re.match('^Sequences producing significant alignments:', line) or re.match('^>', line)):
			if(reading) :
				break
			reading = True
			out = "";
		else:
			out += line

	for line in out.split("\n"):
		if (line != ""):
			exons.append(line.split(" ")[0])
	return exons
	

def parse_descr_file (descr_file) :

	species_name="";
	short_name="";
	params=[];
	all_species = {};
	new_species = 1;
	
	
	FILE = open(descr_file,'r');
	while FILE :
		species_name = FILE.readline().rstrip('\n');
		if (not species_name):	# end of file
			break
		short_name = FILE.readline().rstrip('\n');
		params_line = FILE.readline().rstrip('\n').split(' ')[1];
		params = params_line.split(':');
		FILE.readline()
		all_species[short_name] = (species_name, params);
		
		
	FILE.close()
	return all_species
	
	
# functions to aid writing and formatting the log file
def write_species_info (logfile, short_name, (species_name, params)):
	logfile.write("--SPECIES INFO: " + species_name + " --\n\n");
	logfile.write("Source .fasta file: " + regions_path + short_name + ".fasta\n")
	logfile.write("Type of sequence:   " + params[0] + "\n")
	logfile.write("Sequence ID:        " + params[2] + "\n")
	logfile.write("Region of interest: " + params[3] + " - " + params[4] + "\n")
	logfile.write("Strand:             " + params[5] + "\n\n")

def annotate_exons (exons, short_name):
	outfile = './blastout/'+short_name+'.blastout';
	
	exons_found = {}
	i = 0
	found_f = False
	results = open(outfile, "r")
	for line in results.readlines():
		if(i == len(exons)):break;
		if(re.match('^>{0} '.format(exons[i]), line)):
			found_f = True
		else:
			if(found_f == True):
				reg = re.match(r'^Query: (\d*)', line)
				if(reg):
					exons_found[exons[i]] = reg.group(1)
					found_f = False
					i += 1
	### Anotacija i analiza za neprosirene
	fasta_file = regions_path + short_name + '.fasta'
	file = open (fasta_file, 'r');
	file.readline();
	seq = file.read();
	seq = "".join(seq.split('\n'))

	out = ""
	for iter in range(1, number_of_exons + 1):
		if exons_found.has_key('{0}'.format(iter)):
			out += "  T"
		else:
			#pronadji lijevu granicu
			left = 0
			reversed_range = range(0, iter+1)
			reversed_range.reverse() 
			for iterL in reversed_range:
				#print "iterL:" + str(iterL)
				if exons_found.has_key('{0}'.format(iterL)):
					left = int(exons_found['{0}'.format(iterL)])
					break;
			#pronadji desnu granicu
			right = 0
			for iterR in range(iter, number_of_exons + 2):
				if exons_found.has_key('{0}'.format(iterR)):
					right = int(exons_found['{0}'.format(iterR)])
					break;
			#print "Left:{0}; Right:{1}".format(left, right)
			if(left != 0 and right != 0):
				matchObj = re.findall(r'N+', seq[left:right])
				if len(matchObj) > 0 :
					out += "  N"
				else:
					out += "  ?"
			else:
				if(left != 0):
					matchObj = re.findall(r'N+', seq[left:])
					if len(matchObj) > 0 :
						out += " NX"    
					else:
						out += "  X"
				else: 
					if(right != 0):
						matchObj = re.findall(r'N+', seq[:right])
						if len(matchObj) > 0 :
							out += " NX"    
						else:
							out += "  X"
	return out
	
def write_exon_info (logfile, exons, exons_widened, exons_proteins, short_name) :
	
	if (exons == -1) :
		logfile.write ("No appropriate source file found (was not able to fetch the region from ensembl exon_db).\n")
		return
	####
	#MARIO
	out = ""
	for iter in range(1, number_of_exons + 1):
		out += "{0:3d}".format(iter)
	logfile.write(out)
	logfile.write('\n')
	
	annotation = annotate_exons(exons, short_name)
	logfile.write(annotation)
	logfile.write('\n');
	###
	#statistics_original
	if annotation != "":
		statistics[5] += 1
		all_exon_info = re.split(' *', annotation)
		exon_counter = -2; #prvi u all_exon_info je " "
		for s in all_exon_info:
			exon_counter += 1;
			if s == 'T':
				exon_statistics[exon_counter] += 1;
				statistics[0] += 1
			elif s == 'N':
				statistics[1] += 1
			elif s == 'X':
				statistics[2] += 1
			elif s == 'NX':
				statistics[3] += 1
			elif s == '?':
				statistics[4] += 1
	
		
	# if widened search is provided
	if (exons_widened != -1) :
		annotation = annotate_exons(exons_widened, short_name + "_whole")
		logfile.write(annotation)
		logfile.write('\n');
		###
		#statistics_original
		if annotation != "":
			statistics_w[5] += 1
			all_exon_info = re.split(' *', annotation)
			for s in all_exon_info:
				if s == 'T':
					statistics_w[0] += 1
				elif s == 'N':
					statistics_w[1] += 1
				elif s == 'X':
					statistics_w[2] += 1
				elif s == 'NX':
					statistics_w[3] += 1
				elif s == '?':
					statistics_w[4] += 1
	
	# Biljezenje exona pronadjenih preko proteina
	new_protein_exons = []
	for i in range(1, number_of_exons + 1):
		if str(i) in exons_proteins:
			logfile.write("  T")
		else:
			logfile.write("   ")
		if (str(i) not in exons_proteins) & (str(i) in exons_widened):
			protein_exon_discovery_statistics[i-1] += 1
			new_protein_exons.append(i)
		protein_exon_discovery_species[short_name] = new_protein_exons
	
	
	#/MARIO
	####
	logfile.write ('\n\n\n')
	
	
	
def write_N_region_info(logfile, lenght_of_Ns) :
	if (length_of_Ns != -1) :
		logfile.write ('N Region Status: ');
		if (length_of_Ns == 0) :
			logfile.write ('No N regions found in sequence.\n\n')
		else :
			logfile.write ('Found ' + str(length_of_Ns) + " N regions.\n\n")
		
def new_exon_discovery(species, exons, exons_w):
	exon_dict = {}
	new_exons = []
	for e in exons:
		exon_dict[e] = 1
	for e in exons_w:
		if exon_dict.has_key(e) == False:
			new_exons.append(e)
			new_exon_discovery_statistics[int(e) - 1] += 1
	new_exon_discovery_species[species] = new_exons
			
		
# execute

log = open(log_file, 'w')	
all_species = parse_descr_file(descr_file)
exons={}
exons_w={}
exons_p = {}
for short_name in all_species :

	write_species_info(log, short_name, all_species[short_name]);	# log
	
	(species_name, params) = all_species[short_name];
	
	length_of_Ns = test_for_N_regions(regions_path + short_name + '.fasta');
	write_N_region_info(log, length_of_Ns);							#log
	
	exons[short_name] = find_exons(params, species_name, short_name)
	exons_w[short_name] = find_exons(params, species_name, short_name + '_whole')
	
	exons_p[short_name] = find_exons_prot(params, species_name, short_name)
	if exons_w[short_name] == -1: continue
	
	write_exon_info(log, exons[short_name], exons_w[short_name], exons_p[short_name], short_name);	# log
	
	new_exon_discovery(short_name, exons[short_name], exons_w[short_name])
	
	print short_name
	print "Exons (predefined):  ", exons[short_name]
	print "Exons (widened area):", exons_w[short_name]
	print '\n\n'

log.write("Statistika:\n")	
log.write("izvorno\n")
for i in range(0, 5):
	statistics[i] = float(statistics[i]) / (statistics[5] * number_of_exons) * 100
log.write("T = {0:.2f}% \nN = {1:.2f}% \nX = {2:.2f}% \nNX = {3:.2f}% \n? = {4:.2f}%\n".format(statistics[0], statistics[1], statistics[2], statistics[3], statistics[4]))

log.write("prosireno\n")
for i in range(0, 5):
	statistics_w[i] = float(statistics_w[i]) / (statistics_w[5] * number_of_exons) * 100
log.write("T = {0:.2f}% \nN = {1:.2f}% \nX = {2:.2f}% \nNX = {3:.2f}% \n? = {4:.2f}%\n".format(statistics_w[0], statistics_w[1], statistics_w[2], statistics_w[3], statistics[4]))

log.write("Prosjecni broj pronadjenih exona (%):\n")
out = ""
for iter in range(1, number_of_exons + 1):
	out += "{0:6d}".format(iter)
log.write(out)
log.write('\n')
for i in range(0, number_of_exons):
	exon_statistics[i] = float(exon_statistics[i]) / statistics[5] * 100
	log.write(" {0:05.2f}".format(exon_statistics[i]))
log.write('\n')

log.write("Prosjecni broj exona otkrivenih prosirenjem:\n")
log.write(out)
log.write('\n')
for i in range(0, number_of_exons):
	new_exon_discovery_statistics[i] = new_exon_discovery_statistics[i]/((100-exon_statistics[i])/100*statistics[5]) * 100
	log.write(" {0:05.2f}".format(new_exon_discovery_statistics[i]))
log.write('\n')
log.write('Sto je prosirenje naslo:\n')
for s in new_exon_discovery_species:
	if new_exon_discovery_species[s]:
		log.write("{0} : {1}\n".format(s, new_exon_discovery_species[s]))
log.write('\n')	

log.write('Sto se preko DNA naslo a nije preko proteina:\n')
log.write(out)
log.write('\n')	
for i in range(0, number_of_exons):
	log.write(" {0:5d}".format(protein_exon_discovery_statistics[i]))
log.write('\n')
for s in protein_exon_discovery_species:
	if protein_exon_discovery_species[s]:
		log.write("{0} : {1}\n".format(s, protein_exon_discovery_species[s]))
log.close()
