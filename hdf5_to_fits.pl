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
$filename = $ARGV[0];

&read_flash("$filename");

# Define square matrix
$ncart = 800;
$plot_var_cart = zeroes($ncart,$ncart);

# Add zero values between 0:$xrange[0]
my $nx2 = int ($nx * $xrange[1]/($xrange[1]-$xrange[0]));
my $nxi = $nx2 - $nx;
my $nxf = $nx2 - 1;

$plot_var_new = zeroes("$nx2","$ny");

# zero in 0:$nxi-1
$plot_var_new->slice("$nxi:$nxf,:") .= $plot_var;

# Convert to Cartesian coords
$ts = t_linear(scale => [2.0*pi/($ny-1), 140./$ncart]); #scale in y and x
$tu = !t_radial();

$plot_var_cart = transpose($plot_var_new)->map($tu x $ts,$plot_var_cart,{method=>'l'});
# $plot_var_cart = $plot_var_new->map($tu x $ts,$plot_var_cart,{method=>'l'});

# Transpose and rotate 90 degrees
wfits transpose($plot_var_cart)->map(t_linear(rot=>-90)), "$filename.fits";

# wfits $plot_var_cart, "$filename.fits";
# $plot_var->wfits("$filename.fits");
