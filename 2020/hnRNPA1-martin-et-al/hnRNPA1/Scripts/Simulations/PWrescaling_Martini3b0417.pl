#!/usr/bin/perl -w
# 2019.03.19
# Yong Wang
# adjust the VDW interactions between Protein and Water in Martini3.0.4.1.7
use strict;
use Getopt::Std;

if ($#ARGV == -1) {
    &usage() and exit -1;
}

my $topfile;
my $top2file;
my $theta;

&getopts('I:T:O:', \my %opts);

if (defined $opts{'I'}) {
    $topfile = $opts{'I'};
} else {
    &usage() and exit;
}

if (defined $opts{'T'}) {
    $theta = $opts{'T'};
} else {
    &usage() and exit;
}

if (defined $opts{'O'}) {
    $top2file = $opts{'O'};
} else {
    &usage() and exit;
}

open(TOP, "$topfile") or die "Unable to open $topfile.";
open(TOPNEW, ">$top2file") or die "Unable to open $top2file.";

#=============================
#=============================

while (my $line = <TOP>) {
        if ($line =~ /\[ nonbond_params \]/) {
        printf (TOPNEW "\n");
        printf (TOPNEW "$line");
        last;
        }
        printf (TOPNEW "$line");
}


#=============================
#[ nonbond_params ]
#=============================
while (my $line = <TOP>) {
  if ($line =~ /\[ moleculetype \]/) {
        printf (TOPNEW "\n");
        printf (TOPNEW "$line");
	last;
  }

  if ($line=~/W/) {
     	(my $tmp1, my $at1, my $at2, my $tmp, my $sigma, my $epsilon)=split(/\s+/,$line);
	# Ion bead: TQ5
	# Water beads: W, SW and TW
	if ( $tmp =~ /1/ && ($at1 !~ /W/ || $at2 !~ /W/) && $line !~ /TQ5/) {
     	printf(TOPNEW " %5s %5s  1  %7.6E    %7.6E ; original epsilon %7.3f * %7.3f \n",$at1,$at2,$sigma,$epsilon*$theta,$epsilon,$theta);
	} else {
     	printf(TOPNEW "$line");
	}
  } else {
     	printf(TOPNEW "$line");
  }
}

while (my $line = <TOP>) {
        printf (TOPNEW "$line");
}

close(TOP);
close(TOPNEW);
exit;
#======================================

sub usage {
    print STDERR << "EOF";


EOF
}
