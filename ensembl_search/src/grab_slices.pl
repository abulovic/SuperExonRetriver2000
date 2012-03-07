#!/usr/bin/perl

# očekuju se dva argumenta komandne linije:
# 1. ime datoteke u kojoj se nalazi popis potrebnih regija
# 2. put do direktorija u koji je naknadno potrebno spremiti sve datoteke

use Bio::EnsEMBL::Registry;

if (@ARGV < 2) {
	print "GRESKA!!\nNedovoljan broj ulaznih argumenata.\nPotrebno je unijeti put do datoteke s popisom vrsta i put do direktorija u kom ce se spremati fasta datoteke.\n";
	print "Opcionalni argument: -w (ucitavaju se i zapisuju u fasta datoteke i cijeli scaffoldi/contizi.)";
	exit(-1);
}

# txt dokument s popisom vrsta i potrebnih regija, predan kao komandnolinijski argument
my $file = $ARGV[0];

# spajanje na EnsEMBL bazu
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

# otvaranje datoteke i parsiranje
open FILE, $file or die $!;
my @lines = <FILE>;
my $new_species = 0;
my $slice_adaptor;
my $species_short_name;
my $short_name_file = "./Regions/species_list.txt";
open (SHORT_NAMES_FILE, ">${short_name_file}");
#=comment
for my $line (@lines) {
	chomp($line);
	if (!$line) { # prazan string
		$new_species = 0;
		next;
	}
	
	if ($new_species == 0) {
		@name = split (/_/, $line);
		$animal_name = "${name[0]} ${name[1]}";
		print $animal_name, "\n";
		print $line, "  Getting slice adaptor...\n";
		$slice_adaptor = $registry->get_adaptor($animal_name, 'Core', 'Slice');
		$new_species++;
		next;
	}
	if ($new_species == 1) {
		print $line, "  Remebering short name....\n";
		$species_short_name = $line;
		$new_species++;
		next;
	}
	if ($new_species == 2) {
		print $line, "   params....\n";
		my @tmp = split (/ /, $line);
		my @params = split(/:/, $tmp[1]);
		my $slice =  $slice_adaptor->fetch_by_region ($params[0], $params[2], $params[3], $params[4], $params[5]);
		my $whole_slice = $slice_adaptor->fetch_by_region($params[0], $params[2]);

		if (!$slice){
			print "Slice not found for ", $species_short_name, "\n";
			next;
		}

		# if slice is found, write species short name to species list
		print SHORT_NAMES_FILE "${species_short_name}\n";

		if ($slice->start() > $whole_slice->end()) {
			print $params[0], " ends before the beggining of the requested region.\n";
			next;
		}

		if ($slice->end() > $whole_slice->end()) {
			print $params[0], " shorter then requested region - cutting down from ", $slice->end(), " to ", $whole_slice->end(), "\n";
			$slice = $slice_adaptor->fetch_by_region ($params[0], $params[2], $params[3], $whole_slice->end(), $params[5]);
		}

		if (@ARGV == 3 and $ARGV[2] == "-w") {
			my $fasta_file_name_whole = "${ARGV[1]}/${species_short_name}_whole.fasta";
			open (FASTA_FILE, ">${fasta_file_name_whole}");
			print $fasta_file_name_whole, "\n";
			print FASTA_FILE ">${species_short_name}\n";
			if ($params[0] == "chromosome") {
				$slice_ex = $slice_adaptor->fetch_by_region ($params[0], $params[2], $params[3]-100000, $params[4]+100000, $params[5]);
			}
			else{
				$slice_ex = $whole_slice;
			}
			my $seq = $slice_ex->seq();
			$seq =~ s/(.{1,60})/$1\n/gs;
			print FASTA_FILE $seq;
			close (FASTA_FILE);
		}

		my $fasta_file_name = "${ARGV[1]}/${species_short_name}.fasta";
		open (FASTA_FILE, ">${fasta_file_name}");
		print $fasta_file_name, "\n";
		print FASTA_FILE ">${species_short_name}\n";
		my $seq = $slice->seq();
		$seq =~ s/(.{1,60})/$1\n/gs;
		print FASTA_FILE $seq;
		close (FASTA_FILE);
	}
	
	
}
close (SHORT_NAMES_FILE);
#=cut

=comment
 my $slice_adaptor = $registry->get_adaptor('Otolemur garnettii', 'Core', 'Slice');

 my $slice = $slice_adaptor->fetch_by_region ('scaffold', 'scaffold_117055', '730', '1318', '-1');

if (!$slice) {
	print "Slice not found\n";
}
else {
	 my $sequence = $slice->seq();
}
 print $sequence, "\n";
=cut


