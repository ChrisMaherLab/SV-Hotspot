#!/usr/bin/env perl

# Created by HD <hd@vt.edu> Jul. 2011
# Last updated: Dec. 2014, v2: allow specifying ID column
# v2: pipe input

my $U = '

	cat <filter-file> | thisfile.pl <filter-col> <select-file> <select-col> <remove?:0|1, default=0> <regex?:0|1, default=0> > <output-file>
	
	Select rows from select-file based on filter-file
    v2: pipe input

	If remove = 1 => select, 0 - remove
	If regex = 1 => use regular expression in matching
	

';


use strict;

scalar @ARGV >= 3 or die $U;

my ($col2, $iFile, $idColumn, $remove, $regex) = @ARGV;
defined $remove or $remove = 0;
defined $regex or $regex = 0;

my %id2keep = ();
foreach my $l (<STDIN>){
	chomp $l;
	if ($l =~ /(^\s*#|^\s+$)/){next};
	$l =~ s/(^\s+|\s+$)//g;
    my @ws = split "\t", $l;
	my $id = $ws[$col2];
    $id2keep{$id} = 1;
}


open F, $iFile or die "Cannot read $iFile\n";
my @ls = <F>;
close F;
foreach my $l (@ls){
	chomp $l;
	my @ws = split /\t/, $l;
	my $id = $ws[$idColumn];
	if (! $regex){
	    if ($remove and ! exists $id2keep{$id}){
		    print "$l\n";
	    }elsif(! $remove and  exists $id2keep{$id}){
	            print "$l\n";
	    }
	}else{
	    if ($remove and ! idMatched($id)){
		    print "$l\n";
	    }elsif(! $remove and idMatched($id)){
	        print "$l\n";
	    }
	}
}

sub idMatched{
    my $id2 = shift;
    foreach my $id (keys %id2keep){
        if ($id2 =~ /$id/){
            return 1;
        }
    }
    return 0;
}
