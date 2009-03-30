#!/usr/bin/perl -w
# plot n-body data
# usage: nbody.pl <nboy_file>

use PDL;
use Fortran::Namelist;

die "Please give data file\n" if (! $ARGV[0]);
my ($file_string, $file_number) = $ARGV[0] =~ /^(.*)_([0-9]*)/;
my $end_number = $ARGV[1];
my $step = 1;

# open ($output,">$fileout");
for ($j=$file_number;$j<=$end_number;$j=$j+$step) {
   $number = sprintf("%04u", $j);
   $filename = "$file_string"._."$number";
   $nml=Fortran::Namelist->new( file => $filename );
#    @list = list($r_body);
#    print @list . "\n";
}
# close ("$output");
