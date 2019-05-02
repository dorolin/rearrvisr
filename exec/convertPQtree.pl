#!/usr/bin/perl

## ------------------------------------------------------------------------
## read PQ tree from Anges output (do not use PQR tree) and convert it
## to a table-like format with columns $marker, $orientation, $car,
## and at least two additional columns with node $type* and $elem*

## Example: convertPQtree.pl XXX_PQTREE_HEUR
## ------------------------------------------------------------------------

use warnings;
use strict;

die "USAGE: convertPQtree.pl PQTREE\n" unless @ARGV==1;

my $in = shift(@ARGV);

open (IN, $in) or die "Cannot read $in file\n";

my $out = $in;
$out = $out.".txt";

open (OUT, ">", "$out") or die "Cannot open $out file\n";


my $line = '';
my @tree = ();
my $car = '';
my $elem = '';
my $sign = '';
my @nodestack = ();
my @countstack = ();


while ($line = <IN>){
    chomp $line;
    if($line =~ s/^>//){
	print "\nconverting data for $line";
    }
    elsif($line =~ s/^#CAR//){
	$car = $line;
    }
    else{
	@tree = split(/\s/,$line);
	while (@tree > 0){
	    $elem = shift(@tree);
	    if(@countstack>0){ ## already node(s)
		## -> increment node counter
		$countstack[-1]++;
	    }
	    if($elem =~ m/_[P|Q]/){ ## node opens
		$elem =~ s/_//;
		push(@nodestack,$elem);
		push(@countstack,0);
	    }
	    elsif($elem =~ m/[P|Q]_/){ ## node closes
		$elem =~ s/_//;
		die "Something in the hierarchy of CAR$car is wrong\n" if $nodestack[-1] ne $elem;
		$elem = pop(@nodestack);
		$elem = pop(@countstack);
	    }
	    else{ ## marker
		if($elem =~ s/-//){
		    $sign = '-';
		}
		else{
		    $sign = '+';
		}
		print OUT $elem." ".$sign." ".$car;
		for (my $i = 0; $i<scalar(@nodestack); $i++){
		    print OUT " ".$nodestack[$i]." ".$countstack[$i];
		}
		print OUT "\n";
	    }
	}
    }
}

if($car eq ''){
    print "\n... no CARs found ...\n\n";
}
elsif($car == 1){
    print "\n... processed $car CAR ...\n\n";
}
else{
    print "\n... processed $car CARs ...\n\n";
}



close (IN);
close (OUT);

exit;

