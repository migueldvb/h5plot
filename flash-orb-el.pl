#!/usr/bin/perl -w

use PDL;
use PDL::Graphics::PLplot;
use PDL::NiceSlice;
use PDL::AutoLoader;
use PDL::IO::Misc;
use Getopt::Long qw [:config pass_through];
use Math::Trig qw [pi];

# $fileout='test';
# $mstar = 9.1732784568533237e-05;
# $mplanets = pdl (1.384419495666871e-08, 4.4922999961435205e-09, 7.430659742048716e-09);
# $mu = $mplanets + $mstar;

GetOptions ( "p=i"=> \@planet);

$file_string=$ARGV[0];
if (($ARGV[1]) && ($ARGV[1]>=0) ) {
   $end_number = $ARGV[1];
} else {
   $end_number = 9999;
}

$tend=100;
$step=2;
$dist=0.041;
@planet=(1,2);

cmap1_init();
plinit;
# plParseOpts(\@ARGV,PL_PARSE_SKIP | PL_PARSE_NOPROGRAM);
plenv(0,$tend,0.0002,0.01,0,0) if (@planet == 0);
# plenv(0,$tend,0.9*$dist,2.*$dist,0,0) if (@planet >= 1);
plenv(0,$tend,0.,0.2,0,0) if (@planet >= 1);

for ($j=1;$j<=$end_number;$j=$j+$step) {
   $number = sprintf("%04u", $j);
   $filename = "$file_string"."$number";
   open ($input,"<$filename") || last;
   close ("$input");
   read_nbody($filename);
   $mu = $m_body + $m_body->(0);
#    print "$x,$y\n";
#    print "$vx,$vy\n";
#    pllab("Time [orbits]","a [AU]","");
   pllab("Time [orbits]","e","");
   foreach $i (@planet) {
       $semia = 1./(2./$r_body->($i) - $v2_body->($i)/$mu->($i));
       $planet=$i;
       $ecc_arr = ($v(:,$planet)**2)->sumover*$r(:,$planet)/$mu($planet) -
          (sum($v(:,$planet)*$r(:,$planet)))/$mu($planet) * $v(:,$planet) -
          $r(:,$planet)/($r_body->($planet));
       $ecc = sqrt(($ecc_arr**2)->sumover);
       print "$time_nbody,$semia,$ecc\n";
#        plpoin($time_nbody/2/pi,$semia*$dist,2);
       plpoin($time_nbody/2/pi,$ecc,2);
   }
}
plend;

$fileout=shift @ARGV if ($ARGV[0]);

sub semia { # semi-major axis
   $semia = $mu($1)/(2.*$mu($1)/sqrt($x*$x + $y*$y + $z*$z) - ($vx*$vx + $vy*$vy + $vz*$vz));
   return $semia;
}

sub ecc { # eccentricity
   my $planet = @_;
   $ecc_arr = ($v(:,$planet)**2)->sumover*$r(:,$planet)/$mu($planet) -
      (sum($v(:,$planet)*$r(:,$planet)))/$mu($planet) * $v(:,$planet) -
      $r(:,$planet)/($r_body->($planet));
   $ecc = sqrt(($ecc_arr**2)->sumover);
   return $ecc;
}

sub plot_orbparms {
   my ($planet,$key) = @_;
   semia($planet);
   plenv(0,$time->max,$semia->min,$semia->max,0,0);
   pllab("Time [days]","a [AU]", "Gl581$key");
   plline($time,$semia);
   ecc($planet);
   plenv(0,$time->max,$ecc->min,$ecc->max,0,0);
   pllab("Time [days]","e", "Gl581$key");
   plline($time,$ecc);
}
