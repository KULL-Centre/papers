#!/usr/bin/perl -w
# 2012.2.29
# 2014.5.3
use strict;
use Getopt::Std;

if ($#ARGV == -1) {
    &usage() and exit -1;
}

my $framefile;
my $Moilpdbfile;
my $FFpdbfile;

&getopts('G:F:O:', \my %opts);

if (defined $opts{'G'}) {
    $Moilpdbfile = $opts{'G'};
} else {
    &usage() and exit;
}

if (defined $opts{'F'}) {
    $FFpdbfile = $opts{'F'};
} else {
    &usage() and exit;
}

if (defined $opts{'O'}) {
    $framefile = $opts{'O'};
} else {
    &usage() and exit;
}

#===============================
# Step 1: Read PDB of Moil
#===============================
open(MOILPDB, "$Moilpdbfile") or die "Unable to open $Moilpdbfile.";
my $iatom=0;
my @MOILatomID=();
my @MOILNATCoor=();
my @occupacy=();
my @beta=();
while (my $line = <MOILPDB>) {
        if ( substr($line, 0, 4) =~ /ATOM/) {
		$iatom++;
                my $atomnum = substr($line, 7, 4);
                my $atomtype = substr($line, 12, 4);
		$atomtype=~s/(^\s+|\s+$)//g;
                my $restype = substr($line, 17, 3);
                my $chain = substr($line, 21, 1);
                my $resnum = substr($line, 23, 3);
		my $align=substr($line, 55, 5);
		my $bias=substr($line, 61, 5);
#TOM     13  SD  MET     1      43.420   1.590   8.120  1.00  0.00
#TOM     22  CA  GLU A  23       2.379  13.452  -1.548  1.00  0.00

		$MOILatomID[$iatom][0]=$iatom;
		$MOILatomID[$iatom][1]=$atomtype;
		$MOILatomID[$iatom][2]=$resnum;
		$MOILatomID[$iatom][3]=$chain;
		$MOILatomID[$iatom][4]=$restype;
		$occupacy[$iatom]=$align;
		$beta[$iatom]=$bias;

                my $x = substr($line, 30, 8); $x =~ s/\s+//g;
                my $y = substr($line, 38, 8); $y =~ s/\s+//g;
                my $z = substr($line, 46, 8); $z =~ s/\s+//g;
                $MOILNATCoor[$iatom][0]=$x;
                $MOILNATCoor[$iatom][1]=$y;
                $MOILNATCoor[$iatom][2]=$z;
        }
}
my $NATOM=$iatom;
printf("There are $iatom atoms in total.\n");
close(MOILPDB);

#===============================
# Step 2: Read PDB of FF
#	then build the correlation list
#===============================
open(FFPDB, "$FFpdbfile") or die "Unable to open $FFpdbfile.";
open(OUT, ">$framefile") or die "Unable to open $framefile.";
$iatom=0;
my $atomn=0;
while (my $line = <FFPDB>) {
        if ( substr($line, 0, 4) =~ /ATOM/) {
                my $atomnum = substr($line, 7, 4);
                my $atomtype = substr($line, 12, 4); $atomtype=~s/(^\s+|\s+$)//g;
                my $restype = substr($line, 17, 3);
                my $chain = substr($line, 21, 1);
                my $resnum = substr($line, 23, 3);

		#ILE CD1 in MOIL corresponds to ILE CD in FF
                if ( ($restype =~ /ILE/) && ($atomtype =~ /CD/)) {
                        $atomtype="CD1";
#                        printf("ATOM    %3d  %-3s %3s %1s %3d \n",$iatom,$atomtype,$restype,$chain,$resnum);
                }

#ATOM     13  SD  MET     1      43.420   1.590   8.120  1.00  0.00
         	for (my $i = 1; $i <= $NATOM; $i++) {
			if ( ($MOILatomID[$i][1] eq $atomtype) && $MOILatomID[$i][2]==($resnum+1) && ($MOILatomID[$i][4] eq $restype) ) {
				#printf("ATOM    %3d  %-3s %3s %1s %3d ; %6d\n",$atomnum,$atomtype,$restype,$chain,$resnum,$i);
				printf(OUT "ATOM %6d %3s  %3s   %3d     %7.3f %7.3f %7.3f  %4.2f %4.2f \n",$atomnum,$atomtype,$restype,$resnum,$MOILNATCoor[$i][0],$MOILNATCoor[$i][1],$MOILNATCoor[$i][2],$occupacy[$i],$beta[$i]);
			}
		}
        }
}
close(FFPDB);
#exit;

close(OUT);

#======================================
exit;   

sub usage {
    print STDERR << "EOF";
      
usage: perl $0 -G 010.pdb -F T4L_G.pdb -O frame10.pdb
 
EOF
}
