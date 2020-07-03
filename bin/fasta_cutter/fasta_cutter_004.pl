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
#(c) Frederic Pont 2016

use strict ;

########################################
# Read arguments from the command line
########################################

sub ReadArguments {

  my $fasta = ""  ;
  my $cut_size = "" ;
  my $offset = "" ;

  foreach my $a (0..$#ARGV) {

    # help
    if ($ARGV[0] eq "-h") {
      die "Syntax: extract_fasta.pl -f 'fasta_file' -l 'list' \n";
    }
    elsif ($ARGV[0] eq "-help") {
      &PrintHelp;
    }

    # fasta_file
    elsif ($ARGV[$a] eq "-f") {
      $fasta = $ARGV[$a+1];
    }
    # list
    elsif ($ARGV[$a] eq "-c") {
      $cut_size = $ARGV[$a+1];
    }
    elsif ($ARGV[$a] eq "-o") {
      $offset = $ARGV[$a+1];
    }


  }
  print "Fasta file : <$fasta>\n";
  print "Subsequence size: <$cut_size>\n";
  print "Offset size: <$offset>\n";

  return($fasta, $cut_size, $offset ) ;

}

##########################
#     cutt fasta - 1 CPU for last buffer
##########################
sub cut_sg {
  my ($buf_ref, $results_file, $ls , $os) = @_ ;
  my @buf = @$buf_ref ;
  local $/ = "\n>";

  # my $ls = 100 ; # sequence piece length
  # my $os = 75 ; # sequence offset shift

  open(R, ">results/cut_" . $results_file) || die("Couldn't open results file\n");

  my $index = 0 ;
  my @tit = () ;

  foreach my $l (@buf) {
    chomp $l ;
    my ($tit, $seq) = split(/\n/, $l, 2) ;
    if ($tit !~ /^>/) {                                               # add '>' if it is missing at line start
               $tit = ">" . $tit ;
             }

    $seq =~ s/\n//g;                              # remove new lines in sequence
    my $li = length($seq) ;
    #print $seq, " ", $li, "\n" ;

    # total seq <= sequence piece
    if ($li <= $ls){
      #print "$li <= ", length($l), "\n";
      @tit = split(/ /, $tit, 2 ) ;
      my $tit2 = $tit[0] . "_" . $index . " " . $tit[1] ; # add _1, _2 after the identifier
      $index++ ;
      print R $tit2, "\n" ;
      print R $seq, "\n" ;
    }

    my $r = $li/$os ;

    # total seq > sequence piece
    if ($li > $ls){
        for my $i (1 .. int($r)){
            @tit = split(/ /, $tit, 2 ) ;
            my $tit2 = $tit[0] . "_" . $index . " " . $tit[1] ; # add _1, _2 after the identifier
            print R $tit2, "\n" ;
            print R substr($seq, $os*($i-1), $ls), "\n" ;
            $index++ ;
        }
        if ($r > int($r)){ # il reste un morceau de sequence à la fin
            @tit = split(/ /, $tit, 2 ) ;
            my $tit2 = $tit[0] . "_" . $index . " " . $tit[1] ; # add _1, _2 after the identifier
            $index++ ;
            print R $tit2, "\n" ;
            #print  "substr  $os  " , int($r), "\n" ;
            print R substr($seq, -$ls), "\n" ;                # take the end of the sequence
        }
    }
  }
  close R ;
}


##########################
# MAIN
##########################
sub main {

  my ($fasta, $cut_size, $offset  ) = ReadArguments();

  # results_file name
                       # remove extention
  $fasta =~ /fasta\/(.*)$/ ;
  my $fasta_name = $1 ;
  my $results_file = $fasta_name ;
  print "results_file : <$results_file>\n" ;


  # open fasta file
  open(FASTA, $fasta) || die("Couldn't read file $fasta\n");

  local $/ = "\n>";

  #my $n = 10 ;
  my $ct = 0 ;
  my @buf = () ;
  while (my $l = <FASTA>) {                                                   # read by blocks of  10000 FASTA records
      if ($#buf < 10000-1 ){
        push(@buf , $l) ;
      }
      if ($#buf == 10000-1){
        cut_sg( \@buf, $results_file, $cut_size, $offset);
        @buf = () ; #die ;
        $ct = 0 ;
      }
      if (eof(FASTA)){
        cut_sg( \@buf, $results_file, $cut_size, $offset);
        print "end of file : $ct lines\n"
      }
      $ct++ ;
  }

 close FASTA ;


}

my $start = time;
main();
my $duration = time - $start;
print "Execution time: $duration s\n";
