#!/usr/bin/perl -w
# plot FLASH HDF5 data
# usage: ./read_hdf5.pl <hdf5_file>

use PDL;
use PDL::AutoLoader;
use PDL::Transform;
use PDL::IO::HDF5;
use PDL::IO::FITS;
use PDL::Image2D;
use Getopt::Long qw [:config pass_through];
use Text::Wrap;
use Math::Trig qw [pi];

# Parse options from command line 
GetOptions ("save=s" => \$f_name,
           "help"   => \$help);

if ($help) {
   print (<<EOT);
   $0 options:
    --help                Print out this message

EOT
   push (@ARGV, "-h");
}
# print "@ARGV \n";
# unshift (@ARGV, $0);

# get the file string from the command line
die "Please give HDF5 data file\n" if (! $ARGV[0]);

# parse file string and number
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

# Define square matrix
# $plot_var_cart = zeroes("$nx","$nx");
$plot_var_cart = zeroes(800,800);

# Add zero values between 0:$xrange[0]
my $nx2 = int ($nx * $xrange[1]/($xrange[1]-$xrange[0]));
my $nxi = $nx2-$nx;
my $nxf = $nx2-1;

$plot_var_new = zeroes("$nx2","$ny");

$plot_var_new->slice("$nxi:$nxf,:") .= $plot_var;

$ts  = t_linear(s => [2.0*pi/($ny-1), 40]);
$tu  = !t_radial();

$plot_var_cart = transpose($plot_var_new)->map($tu x $ts,$plot_var_cart,{method=>'l'});

wfits $plot_var_cart, "$filename.fits";
# $plot_var->wfits("$filename.fits");

