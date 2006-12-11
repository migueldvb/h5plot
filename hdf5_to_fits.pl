#!/usr/bin/perl -w
# plot FLASH HDF5 data
# usage: ./read_hdf5.pl <hdf5_file>

use PDL;
use PDL::AutoLoader;
use PDL::Transform;
use PDL::IO::HDF5;
use PDL::IO::FITS;
# use Astro::IO::CFITSIO;
use PDL::Image2D;
use Getopt::Long qw [:config pass_through];
use Text::Wrap;
use Math::Trig qw [pi];

# Parse ptions from command line 
my $ns = 32;

GetOptions ("save=s" => \$f_name,
           "help"   => \$help);

my @notes = ("Set ns about 16");

if ($help) {
   print (<<EOT);
   $0 options:
    --help                Print out this message

EOT
   print (wrap ('', '', @notes), "\n");
   push (@ARGV, "-h");
}
# print "@ARGV \n";
# unshift (@ARGV, $0);

# get the file string from the command line
die "Please give HDF5 data file\n" if (! $ARGV[0]);

# read file string and number
my ($file_string, $file_number) = $ARGV[0] =~ /^(.*)_([0-9]*)/;
my ($file_base) = $ARGV[0] =~ /^(.*)_hdf.*/;

# read the final file number
if (($ARGV[1]) && ($ARGV[1]>=0) ) {
   $end_number = $ARGV[1];
   # FIXME: accept a file name as $ARGV[1]
#    $end_number = $ARGV[1] =~ /^.*([0-9]{,4})/;
} else {
   $end_number = 9999;
}

die "$end_number should be >= $file_number \n" if ($end_number < $file_number);

$filename = "$file_string"._."$file_number";

&read_flash("$filename");

$plot_var_new = zeroes 4.*($plot_var->dims)[0], 4.*($plot_var->dims)[0];

$plot_var_new = $plot_var->map(t_identity,$plot_var_new,{method=>'s',b=>'e'});

wfits $plot_var_new, "$filename.fits";
# $plot_var->wfits("$filename.fits");

