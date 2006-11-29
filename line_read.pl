#!/usr/bin/perl -w
# Reads line profile data

use PDL;
use PDL::AutoLoader;
use PDL::Graphics::PLplot;
use Getopt::Long qw [:config pass_through];

GetOptions ( "png=s"  => \$png_mode);

$nphi = 12;

# loop over phases
for($i=0; $i<$nphi; $i++) {

      if ($png_mode) {
         # define plot file
         my $fileout="profile_$i${png_mode}050.png";
         plsfnam("$fileout");
         print "Output file: $fileout\n";
         plsdev("png");
      } else {
         plsdev("xwin");
      }

      # read profile ASCII data
      ($v_scaled, $line_prof) = rcols "profile_$i.dat";

      plscol0 (0,255,255,255);
      plscol0 (1,0,0,0);
      plspage (0,0,550,400,0,0);
      plinit;
      cmap1_init();
      plenv($v_scaled->(0), max($v_scaled), 0., max($line_prof), 0, 0);
      pllab ("#fr velocity [km/s]", "#gS#u2#d ", "");
      plline($v_scaled,$line_prof);
      plend;

}
