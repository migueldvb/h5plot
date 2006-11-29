#!/usr/bin/perl -w
# Plot analytical ring profiles with data

use PDL;
use PDL::AutoLoader;
use PDL::IO::HDF5;
use PDL::Graphics::PLplot;
use PDL::Image2D;
use Text::Wrap;
use Math::Trig qw [pi];
use Math::Gsl::Sf qw(:Bessel);
use Getopt::Long qw [:config pass_through];

# ($r,$prof1,$prof2,$prof3) = rcols("besseli.dat", 0,1,2,3, { EXCLUDE => '/^#/' } );

# unshift (@ARGV, $0);

plParseOpts (\@ARGV, PL_PARSE_SKIP );

sub prof {
# calculates analytical ring profile
   $npoints=200;
   $dens = zeroes($npoints);
   $rad = zeroes($npoints);
   my $sf = new Math::Gsl::Sf;
   my $r = new Math::Gsl::Sf::Result;

   my $status = $sf->bessel_Inu_e(0.25,1.,$r);
   $tmp = $r->val;
   print "$tmp \n";

   for ($i=0;$i<$npoints;$i++) {
      my $x = 0.2 + 1.6 * $i/($npoints-1.);
      my $status = $sf->bessel_Inu_e(0.25,2.*$x/$tau,$r);
      $dens->index($i) .= 1/$tau/$x**0.25 * exp(-(1.+$x*$x)/$tau) * $r->val;
      $rad->index($i) .= $x;
   }
}

plscol0 (0,255,255,255);
plscol0 (1,0,0,0);

plinit();

plcol0(1);
$ymax = 2.5;
plenv(0.2,1.8,0.,$ymax,0,1);
pllab ("#frr", "#gS", "Viscous ring (#gn = 4.77 #[Ox00D7] 10#u-5#d)");

$jpt = 2;
# @files=(0,10,40,149);
@files=(0,20,50,150,320);
# for($j=0; $j<=26; $j=$j+10) {
foreach $j (@files) {
   $number = sprintf("%03d", $j);
# read file string and number
   my ($file_string, $file_number) = $ARGV[0] =~ /^(.*)_([0-9]*)/;
#    $filename = "/home/miguel/runs/mako/cmp_problem/A_hdf5_chk_0$number";
   $filename = "$file_string"._."0$number";
   &read_flash($filename);
#    $tau = 0.016+$time*12.*$tvisc;
   $tau = $time*12.*$tvisc;
#    $tau = sprintf("%.5f", $tau);
   print "0.016 + $time * 12. * $tvisc = $tau \n";
   plpoin(1.36,$ymax*0.8-$jpt*0.04*$ymax,$jpt);
   plptex(1.4,$ymax*0.8-$jpt*0.04*$ymax,0,0,0,"#gt = ".sprintf("%.3f", $tau));
   plpoin($x,$profile,$jpt);
#    $filename = "/home/miguel/runs/mako/visc_obj/$filename";
#    &read_flash($filename);
#    $jpt++;
#    plpoin($x,$profile,$jpt);
   &prof();
   plline($rad,$dens);
   $jpt++;
}

plend();
