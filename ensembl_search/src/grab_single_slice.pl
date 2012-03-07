#!/usr/bin/perl

# exit statuses:
#	-1: slice not found
#	 1: requested coordinates outside the range of scaffold / contig

use Bio::EnsEMBL::Registry;

if (@ARGV < 3) {
	print "GRESKA!!\nNedovoljan broj ulaznih argumenata.\nPotrebno je unijeti:\n\t1. Puno ime vrste\n\t2. Ime .fasta datoteke\n\t3. Koordinate regije koju je potrebno dohvatiti\n\t4. Opcionalno -w kao naznaku da se u .fasta datoteku upise citav scaffold/contig.";
	exit(-1);
}

# spajanje na EnsEMBL bazu
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

@name = split (/_/, $ARGV[0]);
$animal_name = "${name[0]} ${name[1]}";
print $animal_name, "\n";
print "  Getting slice adaptor for species ", $animal_name, "...\n";

# dohvacanje adaptora
my $slice_adaptor = $registry->get_adaptor($animal_name, 'Core', 'Slice');

#izgled parametara:
#genescaffold:HEDGEHOG:GeneScaffold_8306:15724:27470:1

my @params = split(/:/, $ARGV[2]);
#fetching region of interest1
my $slice =  $slice_adaptor->fetch_by_region ($params[0], $params[2], $params[3], $params[4], $params[5]);
#fetching whole transcript
my $whole_slice = $slice_adaptor->fetch_by_region($params[0], $params[2]);

if (!$slice){
	print "Slice not found for ", $animal_name, "\n";
	exit(-1);
}


print "Slice end: ", $slice->end(), "\n";
print "Whole slice end: ", $whole_slice->end(), "\n";

if ($slice->start() > $whole_slice->end()) {
	print $params[0], " ends before the beggining of the requested region.\n";
	exit(1);
}

if ($slice->end() > $whole_slice->end()) {
	print $params[0], " shorter then requested region - cutting down from ", $slice->end(), " to ", $whole_slice->end(), "\n";
	$slice = $slice_adaptor->fetch_by_region ($params[0], $params[2], $params[3], $whole_slice->end(), $params[5]);
}

if ($ARGV[3] == "-w"){
	$slice = $whole_slice;
}

my $fasta_file_name = "./Regions/${ARGV[1]}.fasta";
open (FASTA_FILE, ">${fasta_file_name}");
print "Saving sequence to file: ", $fasta_file_name, "\n";
print FASTA_FILE ">${ARGV[1]}\n";
my $seq = $slice->seq();
$seq =~ s/(.{1,60})/$1\n/gs;
print FASTA_FILE $seq;
close (FASTA_FILE);


