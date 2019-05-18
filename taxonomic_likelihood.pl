
###################################################################################################################################
#
#
#	taxonomic_likelihood.pl
#
#    	Copyright (C) 2019 Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
# ##################################################################################################################################
# 
# 
# citation:
# 	Chesters. The phylogeny of insects in the data driven era.
# 
# 
# 	WARNING: alpha release, barely working. set up just for one rank, family.
# 
# 
# 
# Input is 1) phylogeny, 2) taxonomic table of members of the phylogeny
#
# Taxonomic table looks like this:
#
# member	parvorder	tribe	subtribe	subfamily	cohort	class	infraorder	subclass	superorder	family	subgenus	genus	order	suborder	superfamily	infraclass
# Abidama	NA	NA	NA	NA	Paraneoptera	Insecta	Cicadomorpha	Pterygota	NA	Cercopidae	NA	Abidama	Hemiptera	Auchenorrhyncha	Cercopoidea	Neoptera
# Abisara	Heteroneura	NA	NA	Hamearinae	Holometabola	Insecta	Neolepidoptera	Pterygota	NA	Riodinidae	NA	Abisara	Lepidoptera	Glossata	Papilionoidea	Neoptera
# Abispa	NA	NA	NA	Eumeninae	Holometabola	Insecta	Aculeata	Pterygota	NA	Vespidae	NA	Abispa	Hymenoptera	Apocrita	Vespoidea	Neoptera
# 
# 
# You can produce this table yourself, or the following Perl script will produce it, 
# using the NCBI taxonomic database.
# Obviously this will only work if the species you are using, are contained in that database.
# 
# perl taxonomic_report.pl -input itol2_tree.pruned -output itol2_tree.pruned.tax_report -node 6960 -newick
# 
# This script then codes the taxa in your tree as characters for input into Raxml
# perl taxonomic_likelihood.pl -newick itol2_tree.pruned -taxa tabulated_taxa -output tax_on_phylo.OUT -use_ranks family genus
# 
# Then calculate taxonomic lnL with Raxml
# raxmlHPC -f e -s tax_on_phylo.RESULT -m BINGAMMA -n testRaxml -t unrooted.nwk
# 
# 
# 
# 
# #############################################################################################################################
# 
# 
# 
# 	CHANGE LOG
# 
# 2019-05-06: Script started. For revising a paper, need some measures of tree quality to compare two trees.
# 		Was in volved in taxonomic congreunce some time ago, might resurect this.
# 		Will try implement a new version using RaxML.
# 
# 
# 
# 
# 
# 
###############################################################################################################################

my $arg_string  = join ' ', @ARGV;



#####################################
read_command_arguments($arg_string);#
#####################################



#######################
read_tree_terminals();#
#######################



##################
read_tax_table();#
##################


 open(LOG , ">tax_on_phylo.LOG");
 open(RESULTS , ">tax_on_phylo.RESULT");


my @which_ranks_to_test = split /\s+/, $ranks_to_test;

my $logprint = "name\t";
for $i(0 .. $#terminal_labels)
	{
	my $tip = $terminal_labels[$i];
	$logprint .= "$coded\t";
	};
# print LOG "$logprint\n";

my $count_chars =0;

foreach my $rank(@which_ranks_to_test)
	{

	my $names_list_at_rank = $names_for_each_rank{$rank};$names_list_at_rank =~ s/\t+$//;$names_list_at_rank =~ s/^\t+//;
	my @names_at_rank = split /\t+/ , $names_list_at_rank;@names_at_rank = sort @names_at_rank;
	$count_chars += (scalar @names_at_rank) + 1;

	my $count_species = scalar @terminal_labels;

	my $print_header = join '	' , @names_at_rank;
	print LOG "col1\t$print_header\n";
	print "\nprocessing data at rank $rank, there are $#names_at_rank names at this rank:\n"; # $names_list_at_rank\n";
	print RESULTS "$count_species $count_chars\n";

	my %print_out_results;my %print_out_results2;

	foreach my $currentname (@names_at_rank) 
		{
		print "\t$currentname\n";
		my $logprint = "$currentname\t";

		##########################################################################
		for $i(0 .. $#terminal_labels)
			{
			my $tip = $terminal_labels[$i];
			my $taxaname = $tax_name_of_each_rank_for_species{$tip}{$rank};
			my $coded = "?";
			if($taxaname =~ /\w/)
				{
				if($taxaname eq "NA")
					{
					$coded = "?";
					}else{
					if($taxaname eq $currentname){$coded = 1}else{$coded = 0};
					};
				}else{
				$coded = "?";
				};
		#	$logprint .= "$coded\t";
			$print_out_results{$tip} .= "$coded\t";$print_out_results2{$tip} .= "$coded";
			};
		##########################################################################


		};

	my @list_keys = keys %print_out_results;@list_keys = sort @list_keys;
	foreach my $key(@list_keys)
		{
		print LOG "$key\t$print_out_results{$key}\n";
		print RESULTS "$key\t1$print_out_results2{$key}\n";
		};

	};


close LOG;close RESULTS;


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub read_tax_table
{

open(FILE, $taxa_file) || die "\nerror cant open file $taxa_file\n";
my $line_count = 0;
while (my $line = <FILE>)
	{
	$line_count++;
	$line =~ s/\n//;$line =~ s/\r//;
	if($line_count == 1)
		{
		print "taxa file header:$line\n";
	#	$line =~ s/^\w+\t//;
				
		@all_ranks = split /\t/, $line; $ranks_present = scalar @all_ranks; print "there are $#all_ranks ranks in your taxon file\n"; 

		}else{

		@taxa_of_current_member = split /\t/, $line; 
		my $current_member = $taxa_of_current_member[0];
		for my $i(1 .. $#taxa_of_current_member)
			{
			my $which_rank = $all_ranks[$i];
			my $tax_name = $taxa_of_current_member[$i]; #	print "$current_member $i $which_rank $tax_name\n";

			unless($tax_name eq "NA")
				{
				$tax_name_of_each_rank_for_species{$current_member}{$which_rank} = $tax_name;
				# store all names at each rank
				unless($names_for_each_rank{$which_rank} =~ /\t$tax_name\t/){$names_for_each_rank{$which_rank} .= "	$tax_name	"};
				};

			};
		};
	
	};



close FILE;


};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub read_tree_terminals
	{
	open(INTREE, $tree_file) || die "\nerror ... input file you have given ($fas_file) cannot be opened. \n";
	while (my $line = <INTREE>)	{if($line =~ /(.+\(.+)/){$newick_string = $1}};
	close(INTREE);
	$newick_string =~ s/ //g;$tree_length = length($newick_string);
	print "\nnewick string of length $tree_length has been read from file:$tree_file\n\n";

	while($newick_string =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/){$count_scientific_notation_branchlengths_removed++}; 
	# remove regular branchlengths: 0.02048
	while($newick_string =~ s/\:\-*\d+\.\d+//){$count_regular_format_branchlengths_removed++}; 
	# remove 0 length branchlengths
	while($newick_string =~ s/\:\d+//){$count_zero_length_branchlengths_removed++};
	while($newick_string =~ s/(\))\d\.\d+/$1/){$count_proportional_node_support_removed++}; 

#############################################################################################################

print "parsing newick string ....\n";

	while ($newick_string =~ /\(([^\(\)]+)\)/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
		{
		my $node = $1; my $nodeID = "INTERNAL_NODE_$interal_node"; my $length_newick_string = length($newick_string);
		my @child_nodes = split /\,/ , $node;
		if(scalar @child_nodes == 1)# an error, internal node with no split
			{
			if($node =~ /^INTERNAL_NODE_\d+$|INTERNAL_NODE_\d+\:\d+\.\d+$/)
			{
			$newick_string =~ s/\(([^\(\)]+)\)/$1/;
		#	print "why only one child?\n";
			}else{
			print "some other .... length_newick_string:$length_newick_string ID:$nodeID\n";
			if($length_newick_string <= 500){print "\tnode:$node\n";};
			$interal_node++;
			$root_node = $nodeID;

			};

		}
	elsif(scalar @child_nodes == 2)
		{
	#	print "bifurcation $#child_nodes\n";

		$child_counts{$nodeID} = $#child_nodes;

		for $i(0 .. $#child_nodes)
			{
			$nodes{$nodeID}{$i} = $child_nodes[$i];$nodes{$child_nodes[$i]}{parent} = $nodeID;
			unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++}
			};

		$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;

		$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
		$interal_node++;
		}else{
		print "POLYTOMIE, child nodes:$#child_nodes, IN:$interal_node\n";

		while($node =~ s/([^\,]+)\,([^\,]+)/INTERNAL_NODE_$interal_node/)
			{
			my $member1 = $1; my $member2 = $2;			
		#	print "\tmember1:$member1 member2:$member2, repl:INTERNAL_NODE_$interal_node\n";
			$child_counts{$nodeID} = 1;
				$nodes{$nodeID}{0} = $member1;$nodes{$member1}{parent} = $nodeID;
				$nodes{$nodeID}{1} = $member2;$nodes{$member2}{parent} = $nodeID;
			if($node =~ /([^\,]+)\,([^\,]+)/)
				{
				print "\t\tnothing\n";
				}else{
				$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;
				print "\t\tnewick string replacment with INTERNAL_NODE_$interal_node\n";
				};
			$root_node = $nodeID;
			$interal_node++;$nodeID = "INTERNAL_NODE_$interal_node";
			};

			
		};
	}#while ($newick_string =~



@terminal_labels = keys %terminals;@terminal_labels = sort @terminal_labels;
$count_terminal_labels = scalar @terminal_labels;

print "counted terminal in tree:$count_terminal_labels
";



	};

######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


# -newick -taxa -output tax_on_phylo.OUT

sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

# -seqfile -output -node -fasta

if($arguments =~ /-newick\s+(\S+)/)
	{
	$tree_file = $1;
	}else{
	die "\ncommand error, \n"
	}

if($arguments =~ /-taxa\s+(\S+)/)
	{
	$taxa_file = $1;
	}else{
	die "\ncommand error, \n"
	}


if($arguments =~ /-output\s+(\S+)/)
	{
	$output_filename = $1;
	}else{
	$output_filename = "tax_on_phylo.OUTPUT";
	}

if($arguments =~ /-use_ranks\s+([^\-]+)/)
	{$ranks_to_test = $1;
	}else{
	die "\ncommand error\n"
	};



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################









