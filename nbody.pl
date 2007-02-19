#!/usr/bin/perl -w
# plot n-body data
# usage: nbody.pl <hdf5_file>

use PDL;
use PDL::AutoLoader;

die "Please give HDF5 data file\n" if (! $ARGV[0]);
my ($file_string, $file_number) = $ARGV[0] =~ /^(.*)_([0-9]*)/;
my $end_number = $ARGV[1];
my $step = 1;
my $fileout="nbody.txt";

# open ($output,">$fileout");
for ($j=$file_number;$j<=$end_number;$j=$j+$step) {
   $number = sprintf("%04u", $j);
   $filename = "$file_string"._."$number";
   read_nbody($filename);
#    print $output $time_nbody->sclr . $r_body . "\n";
   wcols $time_nbody, ">$fileout";
#    @list = list($r_body);
#    print @list . "\n";
   print $time_nbody . list($r_body) . "\n";
}
# close ("$output");
