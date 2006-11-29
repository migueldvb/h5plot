#!/usr/bin/perl -w

use PDL;
use PDL::AutoLoader;
use PDL::IO::HDF5;
use PDL::Image2D;
use PDL::ImageND;
use PDL::Transform;
use PDL::Graphics::PLplot;
use Math::Trig qw [pi];
use Getopt::Long qw [:config pass_through];

GetOptions ("log"   => \$log_mode,
           "polar"  => \$polar_mode,
           "cart"   => \$cart_mode,
           "png=s"  => \$png_mode,
           "dens"   => \$dens_mode,
           "var=s"  => \$var,
           "vect"   => \$vect_mode,
           "block"  => \$block_mode,
           "prof"   => \$prof_mode,
           "rslice=i"=> \$rslice,
           "pslice=i"=> \$pslice,
           "ns=i"   => \$ns,
           "nbin=i"   => \$nbin,
           "point=i"=> \$point,
           "xzoom=f{2}"=> \@xcoord,
           "yzoom=f{2}"=> \@ycoord,
           "save=s" => \$f_name,
           "help"   => \$help);

# push (@ARGV, "-h");
plParseOpts (\@ARGV, PL_PARSE_SKIP | PL_PARSE_NOPROGRAM );

$filename = $ARGV[0];
my ($file_string, $file_number) = $ARGV[0] =~ /^(.*)_([0-9]*)/;
my ($file_base) = $ARGV[0] =~ /^(.*)_hdf.*/;
$n_body_file = $file_base."_plt_n_body_".$number;

# inclination in radians
$inc = 35.*pi/180.;
$nphi = 12;
$radius = 2.;
$radius2 = $radius*$radius;
# to read the velocity components
$vect_mode = 1;

&read_flash("$filename");
$dens = $plot_var;
# &read_n_body("$n_body_file");
$v_proj = zeroes("$nx","$ny");
# ($i_circ,$j_circ) = whichND($z<$radius2);
# indixes of the inner hole

# define distance to plot line of sight
$ns = 8;
# plot arrow pointing through line of sight
$dist = 3.;

# Maximum projected velocity
$vmax = 1.43; # WARNING: harcoded
$vunit = 193.;
# range of projected velocities
$v_range =  -$vmax + 2.*$vmax*sequence($nbin) / ($nbin-1.) if ($nbin);
$v_scaled =  $vunit*$v_range; # WARNING: hardcoded
# $indx = which ($x < $radius/sqrt(2.));
# $indy = which ($y < $radius/sqrt(2.));
$indxmin = rint $nx/8*(4 - sqrt(2.));
$indxmax = rint $nx/8*(4 + sqrt(2.));
$indymin = rint $ny/8*(4 - sqrt(2.));
$indymax = rint $ny/8*(4 + sqrt(2.));
# print "$indxmin $indxmax\n";
# print "$indymin $indymax\n";

for($i=0; $i<$nphi; $i++) {
   $phi = $i/$nphi*2.*pi;
   $v_proj = $temp_velx*cos($phi)*sin($inc) +
   $temp_vely*sin($phi)*sin($inc);
#    $v_proj = $temp_velx*sin($inc);

   if ($prof_mode) {
      if ($png_mode) {
         my $fileout="profile_$i$png_mode$file_number.png";
         plsfnam("$fileout");
         print "Output file: $fileout\n";
         plsdev("png");
      } else {
         plsdev("xwin");
      }

      #  line profile
      $line_prof =  zeroes($nbin) if ($nbin);
      # Cut the density and projected velocity in the inner hole
#       for($ix=$indxmin; $ix<$indxmax; $ix++) {
#          for($iy=$indxmin; $iy<$indxmax; $iy++) {
#             $ind = max ( which ($v_range < $v_proj->slice("$ix,$iy")));
#             # Intensity is proportional to density squared
#             $line_prof->($ind) += ($dens->slice("$ix,$iy"))**2;
#             wcols $v_scaled, $line_prof, "profile_${i}_b.dat";
#          }
#       }


      # get the index of velocity elements in ascending order
      $isort = qsorti $v_proj->slice("$indxmin:$indxmax,$indymin:$indymax")->flat;
      # sort projected velocity and density for binning
      $v_proj_sort = qsort $v_proj->slice("$indxmin:$indxmax,$indymin:$indymax")->flat;
      $dens_sort = $dens->slice("$indxmin:$indxmax,$indymin:$indymax")->flat->index($isort);
      $dens_sort = $dens_sort**2;
      $ymax = max($dens_sort);
      # write to a file
#       wcols $v_proj_sort, $dens_sort, "profile_${i}_sort.dat";
      # bin the data 
      if ($nbin) {
         $v_proj_bin = zeroes($nbin);
         $v_proj_bin = rebin $v_proj_sort, $v_proj_bin;
         $dens_bin = rebin $dens_sort, $v_proj_bin;
         $ymax = max($dens_bin);
         wcols $v_proj_bin, $dens_bin, "profile_${i}_bin.dat";
      }

#       plscol0 (0,255,255,255);
#       plscol0 (1,0,0,0);
      plspage (0,0,640,500,0,0);
      plinit;
      plenv($v_scaled->(0), $v_scaled->($nbin-1), 0., $ymax, 0, 0);
#          max($dens_sort), 0, 0);
      pllab ("#fr velocity", "#gS#u2#d ", "");
#       plline ($v_scaled, $line_prof);
      if ($nbin) {
         # plot the binned data 
         plpoin ($vunit*$v_proj_bin, $dens_bin, 2);
      } else {
         plpoin ($vunit*$v_proj_sort, $dens_sort, 1);
      }
      plend;
   }
 
   my ($zmin, $zmax) =  f2mnmx ($v_proj);
   print "Extrema of $var:\n";
   print "minimum: $zmin, maximum: $zmax \n";
   ($zmin, $zmax) =  (-$vmax,$vmax);
   my $shedge = $zmin + ($zmax - $zmin) * sequence ($ns) / ($ns-1.);
   my $fill_width = 2;
   my $cont_color = 0;
   my $cont_width = 0;

   if ($png_mode) {
      my $fileout="$var$i$png_mode$file_number.png";
      plsfnam("$fileout");
      print "Output file: $fileout\n";
      plsdev("png");
   } else {
      plsdev("xwin");
   }

   plscol0 (0,255,255,255);
   plscol0 (1,0,0,0);
   plspage (0,0,500,500,0,0);
   plinit;
   cmap1_init();

   pladv(0);
   plvpas(0.11, 0.9, 0.1, 0.9, 1);
#    plvpor(0.2, 0.9, 0.2, 0.9);
   plwind($xrange[0], $xrange[1], $yrange[0], $yrange[1]);
#    plenv($xrange[0], $xrange[1], $yrange[0], $yrange[1], 0, 0);
#    plwind($xmin, $xmax, $ymin, $ymax);
   pllab ("#frx / a", "#fry / a", "#frt = $time");
   plshades ($v_proj, $xrange[0], $xrange[1], $yrange[0], $yrange[1],
      $shedge, $fill_width, $cont_color, $cont_width, 0, 0, 0, 0);
   plcol0(1);
   plbox(0,0,0,0,"btcsn","btcsn");

   #plot line of sight
   $x_los = $dist * cos($phi);
   $y_los = $dist * sin($phi);
   plcol0(1);
   pljoin(0.,0.,$x_los,$y_los);
   # print arrow head
   $x_tmp = 0.9*$dist * cos($phi-pi/32.);
   $y_tmp = 0.9*$dist * sin($phi-pi/32.);
   pljoin($x_tmp,$y_tmp,$x_los,$y_los);
   $x_tmp = 0.9*$dist * cos($phi+pi/32.);
   $y_tmp = 0.9*$dist * sin($phi+pi/32.);
   pljoin($x_tmp,$y_tmp,$x_los,$y_los);
#    $x_vect = zeroes(0,0);
#    $y_vect = zeroes(0,0);
#    $x_vect = 3. * cos($phi);
#    $y_vect = 3. * sin($phi);
#    plvect($x_los,$y_los);

   # plot circles
   $xcirc = -$radius + 2.*$radius*sequence(100)/99.;
   $ycirc = sqrt($radius2-$xcirc*$xcirc);
   plcol0(1);
   plline($xcirc,$ycirc);
   plline($xcirc,-$ycirc);
   # read coordinates of primary and secondary
#    $n_body_file = $file_base."_plt_n_body_".$number;
#    print "$n_body_file \n";
#    &read_n_body("$n_body_file");
   $x_body = pdl(-0.5,0.5);
   $y_body = pdl(0.,0.);
   $r_small = .3;
   $xcirc = -$r_small + 2.*$r_small*sequence(100)/99.;
   $r_small = $r_small*$r_small;
   # plot circle around bodies 
   for ($j=0;$j<2;$j++) {
      $ycirc = $y_body->index($j) + sqrt($r_small-$xcirc*$xcirc);
      plcol0(1);
      plline($xcirc+$x_body->index($j),$ycirc);
      $ycirc = $y_body->index($j) - sqrt($r_small-$xcirc*$xcirc);
      plline($xcirc+$x_body->index($j),$ycirc);
   }
   plend;
}
