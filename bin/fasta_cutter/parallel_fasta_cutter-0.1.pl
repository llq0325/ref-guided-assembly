#!/usr/bin/perl -w

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Written by Frederic PONT.
#(c) Frederic Pont  2016

use strict ;

header() ;

# parse config file
my ($nb_cpu, $command, $cut_size, $offset ) = parse_config_file();

# my $nb_cpu = 8 ;

my $data_folder = "fasta";

opendir(DIR, $data_folder) || die "can't opendir $data_folder: $!";

my @data_dir = readdir(DIR) ;

# remove . and .. in file list
my @list_data_files = map { $_  =~ /^\./ ? () : $_  } @data_dir;

# number of files to process
my $lists_number = $#list_data_files +1 ;

 # if number of list is < nb CPU, then  lower the number of CPU
 $lists_number < $nb_cpu ?  $nb_cpu = $lists_number : 1 ;

# add prefix data/ to the file names
my $list_data_files = join(' '. $data_folder . '/', @list_data_files) ;
$list_data_files = $data_folder . '/' . $list_data_files ;

print "\nFiles to process : $list_data_files\n" ;

my $parallel_cmd = "parallel --bar -j". $nb_cpu ." --no-notice " . $command . " -c $cut_size -o $offset -f ::: $list_data_files";

print "\n$parallel_cmd\n\n" ;

system("$parallel_cmd") ;

header() ;



sub header {
  print "===============================================\n" ;
  print " Fasta Cutter - (c) Frederic PONT 2016\n" ;
  print " Free Software - GNU General Public License\n" ;
  print "===============================================\n" ;

}


#########################
#        parse configuration file
#########################

sub parse_config_file{

	my $config_file = "config.txt" ;
	my $nb_cpu = 2 ;
	my $command = '' ;
	my $key ;

	open(CONF, $config_file);

	while (my $line = <CONF>){
		chomp $line ;
		next if $line =~/^(\s*(#.*)?)?$/;   # skip blank lines and comments;

		if ($line =~/^nb_cpu=/ ){
			($key, $nb_cpu) = split(/=/, $line) ;
			print "calculation on $nb_cpu cpu\n" ;
		}

		if ($line =~/^command=/ ){
			($key, $command) = split(/=/, $line) ;
			print "command :  $command\n" ;
		}

    if ($line =~/^cut_size=/ ){
      ($key, $cut_size) = split(/=/, $line) ;
      print "Subsequence size :  $cut_size\n" ;
    }

    if ($line =~/^offset=/ ){
      ($key, $offset) = split(/=/, $line) ;
      print "Offset :  $offset\n" ;
    }
	}

  if ($offset > $cut_size){
    die "Offset should be <= Subsequence size \n" ;
  }

	return ( $nb_cpu, $command, $cut_size, $offset );
}
