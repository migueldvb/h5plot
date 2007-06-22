#!/usr/bin/perl -w
# plot FLASH HDF5 data
# usage: ./plot_flash.pl <hdf5_file>

use PDL;
use PDL::AutoLoader;
use PDL::Transform;
use PDL::IO::HDF5;
use PDL::Image2D;
use PDL::Graphics::PLplot;
# use PDL::NiceSlice;
use Getopt::Long qw [:config pass_through];
use Text::Wrap;
use Math::Trig qw [pi];

my $ns = 32;
# my $dist = .041;
my $dist = 1.;

# Parse options from command line 
GetOptions ("log"   => \$log_mode,
           "polar"  => \$polar_mode,
           "cart"   => \$cart_mode,
           "png"    => \$png_mode,
           "var=s"  => \$var,
           "vect"   => \$vect_mode,
           "block"  => \$block_mode,
           "stream" => \$stream_mode,
           "prof"   => \$prof_mode,
           "rslice=i"=> \$rslice,
           "pslice=i"=> \$pslice,
           "ns=i"   => \$ns,
           "nrv=i"  => \$nrv,
           "nav=i"  => \$nav,
           "point=i"=> \$point,
           "plane=i{2}"=> \@plane,
           "scale=f{2}"=> \@scale,
           "xzoom=f{2}"=> \@xcoord,
           "yzoom=f{2}"=> \@ycoord,
           "dist=f"    => \$dist,
           "rot=f"  => \$rot,
           "nbody"  => \$nbody,
           "save=s" => \$f_name,
           "help"   => \$help);

my @notes = ("Set ns about 16");

if ($help) {
   print (<<EOT);
   $0 options:
    --log                 Log scale
    --polar               Polar coordinates
    --cart                Cartesian coordinates
    --png                 PNG output
    --var unk             Plot unknown variable
    --vect                Overplot velocity vectors
    --block               Draw blocks
    --ns levels           Set number of shade levels
    --nxv r spacing       Sets spacing of vectors
    --nyv phi spacing     Sets spacing of vectors
    --font number         Selects stroke font set (0 or 1, def:1)
    --help                Print out this message

EOT
   print (wrap ('', '', @notes), "\n");
   push (@ARGV, "-h");
}
# print "@ARGV \n";

# unshift (@ARGV, $0);

# plParseOpts (\@ARGV, PL_PARSE_PARTIAL);
plParseOpts (\@ARGV, PL_PARSE_SKIP | PL_PARSE_NOPROGRAM );

# Calculate streamlines
sub stream($pxinit,$pyinit,$dirt) {
   # Read initial position and direction
   my ($pxinit,$pyinit,$dirt) = @_;
   my $segment_nr = 4000;
   my $iter_nr = 100;
   my $i;
   my $j;
#    my $segment_nr = 40;
#    my $iter_nr = 10;
   my $dt = 0.01;
   my $drint= 0.01;

   $rx = $xrange[0] +($pxinit+0.5)*$dx_fine;
   $ry = $yrange[0] +($pyinit+0.5)*$dy_fine;
   print "rx,ry = $rx $ry \n";

#    $xpos = pdl[$rx,$rx];
#    $ypos = pdl[$ry,$ry];
   @xpos = ($rx,$rx);
   @ypos = ($ry,$ry);

   for($i=1;$i<$segment_nr;$i++) {
      for($j=1;$j<$iter_nr;$j++) {
         get_velocity($rx,$ry,@vel);
#          print "vel = @vel \n";

         if ($geometry eq "polar") {
            $ry = $ry + $dirt*$dt*$vel[1]/$rx;
         } elsif ($geometry eq "Cartesian") {
            $ry = $ry + $dirt*$dt*$vel[1];
         }
         $rx = $rx + $dirt*$dt*$vel[0];

#          print "vel @vel\n";
         last if (( $drint < sqrt(($xpos[1]- $rx)**2+($ypos[1]- $ry)**2) && ($geometry eq "Cartesian")));
         last if (( $drint < sqrt(($xpos[1]- $rx)**2+($xpos[1]*$ypos[1]- $rx*$ry)**2) &&
               ($geometry eq "polar")));
      }

#       print "rx,ry,vel = $rx $ry @vel \n";
      $xpos[0] = $xpos[1];
      $xpos[1] = $rx;
      $ypos[0] = $ypos[1];
      $ypos[1] = $ry;

      last if (($rx < $xrange[0]) || ($rx > $xrange[1]));
      last if ((($ry < $yrange[0]) || ($ry > $yrange[1]-0.5*$dy_fine)) && ($geometry eq "Cartesian"));
#       last if ((($ry < $yrange[0]) || ($ry > $yrange[1]-0.5*$dy_fine)) && ($geometry eq "polar"));
#       last if ((0.5*$drint) > sqrt(($xpos[0]- $rx)**2+($ypos[0]- $ry)**2));
      last if (((0.5*$drint) > sqrt(($xpos[0]- $rx)**2+($xpos[0]*$ypos[0]- $rx*$ry)**2)) && ($geometry eq "polar"));
      if (($ry < $yrange[0]) && ($geometry eq "polar")) {
         $ry = 2.*pi+$ry;
      }
      if (($ry > $yrange[1]) && ($geometry eq "polar")) {
         $ry = $ry-2.*pi;
      }
      @ypos = ($ry,$ry) if (! $cart_mode);
 
#       print "rx @xpos @ypos \n";
      if ($cart_mode) {
         $xcart[0] = $xpos[0]*cos($ypos[0]);
         $ycart[0] = $xpos[0]*sin($ypos[0]);
         $xcart[1] = $xpos[1]*cos($ypos[1]);
         $ycart[1] = $xpos[1]*sin($ypos[1]);
         plline(pdl(@xcart),pdl(@ycart));
      } else {
         plline(pdl(@xpos),pdl(@ypos));
      }
   }
}

# Interpolate the velocity at a given position
sub get_velocity {
   ($rxg,$ryg,@vel) = @_;
   $ix = int(($rxg - $xrange[0] - 0.5*$dx_fine)/$dx_fine);# + 1;
   $qx = ($rxg - $xrange[0] - ($ix+0.5)*$dx_fine)/$dx_fine;

   $iy = int((($ryg - $yrange[0] - 0.5*$dy_fine)/$dy_fine));# + 1;
   $qy = ($ryg - $yrange[0] - ($iy+0.5)*$dy_fine)/$dy_fine;

   # Interpolate to calculate the velocity
   $ix1 = $ix+1;
   $iy1 = $iy+1;
   $iy = $iy%$ny if (($iy >= $ny) && ($geometry eq "polar"));
   $iy = $iy+$ny if (($iy < 0) && ($geometry eq "polar"));
   $iy1 = $iy1%$ny if (($iy1 >= $ny) && ($geometry eq "polar"));
   $iy1 = $iy1+$ny if (($iy1 < 0) && ($geometry eq "polar"));
   $ix = $ix%$nx if (($ix >= $nx) && ($geometry eq "polar"));
   $ix = $ix+$nx if (($ix < 0) && ($geometry eq "polar"));
   $ix1 = $ix1%$nx if (($ix1 >= $nx) && ($geometry eq "polar"));
   $ix1 = $ix1+$nx if (($ix1 < 0) && ($geometry eq "polar"));

#    print "nx,ny $nx $ny \n";
#    print "ix,iy $ix $iy \n";
#    print "ix1,iy1 $ix1 $iy1 \n";
   $vel[0] = $temp_velx($ix,$iy)->sclr * (1-$qx)*(1-$qy) +
        $temp_velx($ix1,$iy)->sclr * $qx*(1-$qy) +
        $temp_velx($ix1,$iy1)->sclr * $qx*$qy +
        $temp_velx($ix,$iy1)->sclr * (1-$qx)*$qy;

   $vel[1] = $temp_vely($ix,$iy)->sclr * (1-$qx)*(1-$qy) +
        $temp_vely($ix1,$iy)->sclr * $qx*(1-$qy) +
        $temp_vely($ix1,$iy1)->sclr * $qx*$qy +
        $temp_vely($ix,$iy1)->sclr * (1-$qx)*$qy;
#    print "vel @vel \n";
}

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

# read the step
if ($ARGV[2]) {
   $step = $ARGV[2];
} else {
   $step = 1;
}

# find initial density using h5dump
my $file0 = "$file_string"._."0000";
my $dens0 = `h5dump -a 'dens/minimum' $file0 | grep "(0)"`;
chomp($dens0);
# keep the digits after the colon 
$dens0 =~ s/^.*: *(\d*)\s*/$1/;
print "dens0 = $dens0 \n";

# find tensor viscosity using h5dump
my $tvisc = `h5dump -d 'real runtime parameters' $file0 | grep -A 1 tvisc | grep "e-"`;
chomp($tvisc);
print "tvisc = $tvisc \n" if ($tvisc);

# plinit;
# cmap1_init();
# loop over files
for ($j=$file_number;$j<=$end_number;$j=$j+$step) {
   $number = sprintf("%04u", $j);
   $filename = "$file_string"._."$number";

   &read_flash("$filename");
   if ($ndim == 3) {
      if ($plane[0] == 0) {
         $plot_var = zeroes("$nx","$ny");
         $plot_var = $plot_var_3d->slice(":,:,$plane[1]");
         print "z = ".$z->index($plane[1])."\n";
      } elsif ($plane[0] == 1) {
         $plot_var = zeroes("$nx","$nz");
         $plot_var = $plot_var_3d->slice(":,$plane[1],:")->clump(2);
         print "y = ".$y->index($plane[1])."\n";
      } elsif ($plane[0] == 2) {
         $plot_var = zeroes("$ny","$nz");
         $plot_var = $plot_var_3d->slice("$plane[1],:,:")->clump(2);
         print "x = ".$x->index($plane[1])."\n";
      }
#       print $z."\n";
   }

# Move the grid in azimuth
   if ($rot) {
      # calculate azimuthal position of planet
      $ypl = ($rot * $time);
      $ypl = $ypl - floor($ypl/2./pi)*2.*pi if ($ypl > 2*pi);
      # Find closest integer
      $pl_ind = rint($ypl/$dy_fine);

      print "ypl = $ypl \n";
      print "pl_ind = $pl_ind $ny/2 \n";
      $plot_var_new = zeroes("$nx","$ny");
      $ny1 = $ny-1;
      if ($pl_ind >= $ny/2) {
         $pl_ind = $pl_ind - $ny/2;
         $tmp1 = $ny/2-$pl_ind;
         $tmp2 = $ny/2+$pl_ind-1;
         print "= $tmp1 $tmp2 \n";
         $plot_var_new->slice(":,0:$tmp1") .= $plot_var->slice(":,$tmp2:$ny1");
         $tmp1 = $tmp1+1;
         $tmp2 = $tmp2-1;
         $plot_var_new->slice(":,$tmp1:$ny1") .= $plot_var->slice(":,0:$tmp2");
      } else {
         $pl_ind = $ny/2 + $pl_ind;
         $tmp1 = 3*$ny/2 - $pl_ind;
         $tmp2 = $pl_ind-$ny/2-1;
         print "= $tmp1 $tmp2 \n";
         $plot_var_new->slice(":,0:$tmp1") .= $plot_var->slice(":,$tmp2:$ny1");
         $tmp1 = $tmp1+1;
         $tmp2 = $tmp2-1;
         $plot_var_new->slice(":,$tmp1:$ny1") .= $plot_var->slice(":,0:$tmp2");
      }
      $plot_var = $plot_var_new;
   }

   $xmin = $xrange[0];
   $xmax = $xrange[1];
   $ymin = $yrange[0];
   $ymax = $yrange[1];
   my $abscissa = 'x';
   my $ord = 'y';
   if ($ndim == 3) {
      if ($plane[0] == 1) {
         $ymin = $zrange[0];
         $ymax = $zrange[1];
         $ny = $nz;
         $abscissa = 'x';
         $ord = 'z';
      } elsif ($plane[0] == 2) {
         $xmin = $yrange[0];
         $xmax = $yrange[1];
         $ymin = $zrange[0];
         $ymax = $zrange[1];
         $nx = $ny;
         $ny = $nz;
         $abscissa = 'y';
         $ord = 'z';
      }
   }

   # calculate boundaries
#    print "@xrange \n";
#    $xmin = $xcoord[0] if (($xcoord[0] > $xrange[0]) && ($xcoord[0] < $xrange[1]));
#    $xmax = $xcoord[1] if (($xcoord[1] < $xrange[1]) && ($xcoord[1] > $xrange[0]));
#    $ymin = $ycoord[0] if (($ycoord[0] > $yrange[0]) && ($ycoord[0] < $yrange[1]));
#    $ymax = $ycoord[1] if (($ycoord[1] < $yrange[1]) && ($ycoord[1] < $yrange[0]));
#    $ymin = $ycoord[0] if ($ycoord[0]);
#    $ymax = $ycoord[1] if ($ycoord[1]);

   # transformation to cartesian
#    $dx_fine = float($xrange[1]-$xrange[0])/$nx;
   $dx_fine = float($xmax - $xmin) / $nx;
#    $dy_fine = float($yrange[1]-$yrange[0])/$ny;
   $dy_fine = float($ymax - $ymin) / $ny;
#    $dy_fine = (2 * pi)/($ny-1.e0);#+5.e-2);
#    $x = $dist * ((xvals("$nx")+.5)*$dx_fine + $xrange[0])->dummy(1,"$ny"+1);
   $x = $dist * ((xvals("$nx")+.5)*$dx_fine + $xmin)->dummy(1,"$ny"+1);
#    $x = ((xvals("$nx")+.5)*$dx_fine + $xrange[0])->dummy(1,"$ny"+1);
#    $y = ((xvals("$ny"+1)+.5)*$dy_fine + $yrange[0])->dummy(0,"$nx");
   $y = ((xvals("$ny"+1)+.5)*$dy_fine + $ymin)->dummy(0,"$nx");

# define a new piddle with one extra row to make transformation
# to Cartesians
   $plot_var_long = zeroes("$nx","$ny"+1);
   $tmp = $ny-1;
   $plot_var_long->slice(":,0:$tmp") .=  $plot_var->slice(':,:');
   $plot_var_long->slice(":,$ny") .=  $plot_var->slice(':,0');
#    print "$yrange[0]\n";
#    print "$yrange[1]\n";
   $cgrid2 = plAlloc2dGrid ($x * cos ($y), $x * sin ($y));

   # log scale
   
   $plot_var->inplace->log10 if ($log_mode);
   $plot_var_long->inplace->log10 if ($log_mode);
#    $plot_var = $plot_var/($dist * ((xvals("$nx")+.5)*$dx_fine + $xmin)->dummy(1,"$ny"))**(-1.5);
#    $plot_var_long = $plot_var_long/$x**(-1.5);


   # PLplot

   if ($png_mode) {
#       my $fileout="$var$png_mode$number.png";
      my $fileout="$var.$number.png";
      $dev = "png";
      plsfnam("$fileout");
      print "Output file: $fileout\n";
      plsdev("$dev");
#    } else {
#       plsdev("xwin");
   }

   print "(n_x,n_y) = ($nx,$ny) \n";
   # calculate time in orbits
   $orbit = $time/2/pi;
   $orbit = sprintf("%.3f", $orbit);
#    $time3f = sprintf("%.3f", $time*2*pi);
   $time3f = sprintf("%.3f", $time);
   print "t = $orbit orbits\n";

   cmap1_init();
# plot averaged profile
   if ($prof_mode) {
# FIXME: Plot several profiles on the same figure
# it should be easy to do
      if ($tvisc) {
# calculate viscous time scale
         $tau = $orbit*12.*$tvisc;
         $tau = sprintf("%.5f", $tau);
         print "tau = $tau \n";
      }
      $y0 = min($profile);
      $y1 = max($profile);
      print "Extrema of averaged $var:\n";
      print "minimum: $y0, maximum: $y1 \n";
      if ($y0 == $y1) {next;}
#       plsdev("null");
      plinit;
#       plenv ($xrange[0], $xrange[1], $y0, $y1, 0, 0);
      plenv ($xmin, $xmax, $y0, $y1, 0, 0);
      if ($tau) {
         pllab ("#frr", "#gS", "#frt = 0.016 + ".$tau." t#dv");
      } else {
         pllab ("#frr", "#gS", "#frt = $orbit orbits");
      }
      if ($point) {
         plpoin($x,$profile,$point);
      } else {
#          plcol0(($j+1)%14);
         plline($x,$profile);
      }
#       plline($x,.0560/$x**2);
      plend;
      next;
   }

   if ($pslice) { # take a slice at a given phi
      $var_slice = $plot_var->slice(":,$pslice");
      $xslice = (xvals("$nx")+.5)*$dx_fine + $xrange[0];
   }

   if ($rslice) {
      # take a radial slice and exchange x and y coordinates
      $var_slice = $plot_var->slice("$rslice,:")->xchg(0,1);
      $xslice = (xvals("$ny")+.5)*$dy_fine + $yrange[0];
   }

   if ($rslice || $pslice) {
      $v_rot = 4.777;
#       $var_slice = $var_slice*$v_rot if ($var eq "velx");
      $y0 = min($var_slice);
      $y1 = max($var_slice);
      print "Extrema of slice of $var:\n";
      print "minimum: $y0, maximum: $y1 \n";
      next if ($y0 == $y1);
      plinit;
      if ($pslice) {
         $tmp = sprintf("%.3f", ($pslice+.5)*$dy_fine/pi);
#          plenv ($xrange[0], $xrange[1], $y0, $y1, 0, 3);
         $xslice = ($xslice - 1.)*70;
         plenv (0., 7., -6.4, -4.2, 0, 0);
#          print "var = $var \n";
#          pllab ("#frr", "$var", "#gf = $tmp #gp");
            pllab ("#frr [AU]", "log(#gS)", "");
#          if ($var eq "dens") {
#             pllab ("#frr", "#gS", "#gf = $tmp #gp");
#          } elsif ($var eq "velx") {
#             pllab ("#frr", "#frv#dr#u [km s#u-1#d]", "#gf = $tmp #gp");
# #             plline($xslice,$v_rot*4.*sqrt(3.)*(log($xslice/$xmin))*(1./3.));
#          } elsif ($var eq "vely") {
#             pllab ("#frr", "#fr#d#gp#u", "#gf = $tmp rad");
#          }
      } elsif ($rslice) {
         $tmp = sprintf("%.3f", ($rslice+.5)*$dx_fine + $xrange[0]);
         plenv ($yrange[0], $yrange[1], $y0, $y1, 0, 0);
         pllab ("#gf", "#fr$var", "#frr = $tmp");
      }
      if ($point) {
         plpoin($xslice,$var_slice,$point);
      } else {
#          plcol0(($j+1)%14);
         plline($xslice,$var_slice);
         plline($xslice,sqrt(1.e-2/$xslice));
#          plline($xslice,$y1*(($xrange[0]+.5*$dx_fine)/$xslice));
      }
      plend;
      next;
   }

   # my $pl = PDL::Graphics::PLplot->new; (DEV => "png", FILE => "dens00$_.png");
   # my $pl = PDL::Graphics::PLplot->new (DEV => "xwin", FILE => ":0");

   my ($zmin, $zmax) =  f2mnmx ($plot_var);
   print "Extrema of $var:\n";
   print "minimum: $zmin, maximum: $zmax \n";

   if (@scale) {
      $indcol = (@scale[1] - $zmin) / ($zmax - $zmin);
      ($zmin, $zmax) =  @scale;
      print "minimum: $zmin, maximum: $zmax \n";
   }

   my $shedge = $zmin + ($zmax - $zmin) * sequence ($ns) / ($ns-1.);

   plscol0 (0,255,255,255);
   plscol0 (1,0,0,0);
   my $ncols = 128;  # set when PALETTE set

   plspage(0,0,800,800,0,0) if ($cart_mode);
   plinit;
   cmap1_init();

   my $fill_width = 2;
   my $cont_color = 0;
   my $cont_width = 0;

   # make contour plot
   pladv(0);
   if ($polar_mode || ! $cart_mode) {
   #    plenv ($xrange[0], $xrange[1], $yrange[0], $yrange[1], 0, 0);
   #    plvsta();
      if ($geometry eq "Cartesian") {
         plvpas(0.2, 0.9, 0.22, 0.9, 1);
   #       plvpor(0.2, 0.9, 0.2, 0.9);
#          plwind($xrange[0], $xrange[1], $yrange[0], $yrange[1]);
         plwind($xmin, $xmax, $ymin, $ymax);
#          plwind(-2,2,-2,2);
#          plwind(-40, 40, -40, 40);
#          pllab ("#fr$abscissa [R#dH#u]", "#fr$ord [R#dH#u]", "#frt = $time3f");
         pllab ("#frx / a", "#fry / a", "#frt = $time3f");
      } else {
         plvpor(0.2, 0.9, 0.2, 0.9);
#          plwind($xrange[0], $xrange[1], $yrange[0], $yrange[1]);
         plwind($xmin, $xmax, $ymin, $ymax);
#          plwind(-2,2,-2,2);
         pllab ("#frr", "#frAzimuth", "#frt = $orbit");
      }
      if ($zmin == $zmax) {
         plcol1(0);
         $xdomain = pdl[$xrange[0],$xrange[0],$xrange[1],$xrange[1]];
#          $y = pdl[$yrange[0],$yrange[1],$yrange[1],$yrange[0]];
         $ydomain = pdl[$ymin,$ymax,$ymax,$ymin];
         plfill ($xdomain,$ydomain);
      } else {
         if (@scale) {
#             plcol1($indcol);
            plcol1(1);
            $xdomain = pdl[$xmin,$xmin,$xmax,$xmax];
            $ydomain = pdl[$ymin,$ymax,$ymax,$ymin];
            plfill ($xdomain,$ydomain);
         }
#          plshades ($plot_var, $xrange[0], $xrange[1], $yrange[0], $yrange[1],
         plshades ($plot_var, $xmin, $xmax, $ymin, $ymax,
            $shedge, $fill_width, $cont_color, $cont_width, 0, 0, 0, 0);
      }
   } elsif ($cart_mode) {
#       plenv (-$xrange[1], $xrange[1], -$xrange[1], $xrange[1], 0, 0);
#       plvpor(0.2, 0.9, 0.1, 0.9);
      ($zmin, $zmax) =  f2mnmx ($plot_var_long);
      ($zmin, $zmax) =  @scale if (@scale);
      $shedge = $zmin + ($zmax - $zmin) * sequence ($ns) / ($ns-1.);
 
      plvpas(0.2, 0.9, 0.2, 0.9, 1);
      if ($xrange[1] <= 0.1*pi/2.) {
         $x_min = $xrange[0]*cos($yrange[1]);
         $y_min = $xrange[1]*cos($yrange[1]);
         print "$x_min $y_min \n";
         plwind($x_min, $xrange[1], -$y_min, $y_min);
#          plwind(-2,2,-2,2);
      } else {
#          plwind(-$xrange[1], $xrange[1], -$xrange[1], $xrange[1]);
         plwind (-$xrange[1]*$dist, $xrange[1]*$dist, -$xrange[1]*$dist, $xrange[1]*$dist)
#          plwind(-2,2,-2,2);
      }
#       plwind($xmin, $xmax, $ymin, $ymax);
      pllab ("#frx/AU", "#fry/AU", "#frt = $orbit");
#       pllab ("#frx", "#fry", "");

# for the first data file with uniform density
      if ($zmin == $zmax) {
#          $x = pdl[-$xrange[1],-$xrange[1],$xrange[1],$xrange[1]];
#          $y = pdl[-$xrange[1],$xrange[1],$xrange[1],-$xrange[1]];
# plot torus between $xrange[0] and $xrange[1] with color pcol1(0)
         $nseq = 50;
         $xseq = $xrange[1] * cos(2. * pi/$nseq*sequence($nseq + 1));
         $yseq = $xrange[1] * sin(2. * pi/$nseq*sequence($nseq + 1));
         plcol1(0);
         plfill ($xseq,$yseq);
         $xseq = $xrange[0] * cos(2. * pi/$nseq*sequence($nseq + 1));
         $yseq = $xrange[0] * sin(2. * pi/$nseq*sequence($nseq + 1));
         plcol0(0);
         plfill ($xseq,$yseq);
      } else {
         plshades ($plot_var_long, -$xrange[1]*$dist, $xrange[1]*$dist, -$xrange[1], $xrange[1],
            $shedge, $fill_width, $cont_color, $cont_width, 0, 0,
            \&pltr2, $cgrid2);
      }
      # Plot eccentric orbit
      $ecc = 0.2;
      $a = 1.*$dist;
#       $ecc = 0.2;
#       $a = 6.098*$dist;
#       $b = $a * sqrt(1.-$ecc**2);
#       $xtmp = -(1+$ecc)*$a + sequence(1001)*$a/500;
#       $ytmp = $b * sqrt(1.-($xtmp+$ecc*$a)**2/$a**2);
      $phitmp = -pi + 2*pi*sequence(1001)/1000;
      $rtmp = $a * (1. - $ecc**2)/(1. + $ecc*cos($phitmp));
      $xtmp = $rtmp * cos($phitmp);
      $ytmp = $rtmp * sin($phitmp);
#       plline($xtmp,$ytmp);
#       plline($xtmp,-$ytmp);
   }

# plot streamlines
   if ($stream_mode) {
      plcol0 (1);
      for ($i=0;$i<$ny;$i=$i+10) {
         $tmp = int(0.9/1.9*$nx);
         &stream($tmp,$i,1);
         &stream($tmp,$i,-1);
#          next if $x($tmp)
#          $tmp = int(1.5/1.9*$nx);
#          &stream($tmp,$i,1);
#          &stream($tmp,$i,-1);
      }
      for ($i=2;$i<$nx;$i=$i+10) {
         &stream($i,$ny/2,1);
         &stream($i,$ny/2,-1);
      }
#             &stream(40,40,1);
   }

   plcol0(1);

   # print title
#    pllab ("", "", "#fr$var");
#    if ($var = "dens") {
#       pllab ("", "", "#frSurface density");
#    } 
   
   plbox(0,0,0,0,"btcsn","btcsn");

   # draw blocks
   if ($block_mode && $cart_mode) {
      for($i=0; $i<$num_plot_blocks; $i++) {
         $cur_blk = $index_good->index($i);
         # find boundaries of current block
         $r0 = $bnd_box->slice("0,0,$cur_blk")->sclr;
         $r1 = $bnd_box->slice("1,0,$cur_blk")->sclr;
         $p0 = $bnd_box->slice("0,1,$cur_blk")->sclr;
         $p1 = $bnd_box->slice("1,1,$cur_blk")->sclr;
         # plot radial spokes
         pljoin($r0*cos($p1),$r0*sin($p1),$r1*cos($p1),$r1*sin($p1));
         pljoin($r0*cos($p0),$r0*sin($p0),$r1*cos($p0),$r1*sin($p0));
         # plot azimuthal arcs
         $xseq = cos($p0+($p1-$p0)/50*sequence(51));
         $yseq = sin($p0+($p1-$p0)/50*sequence(51));
         $xin = pdl($r0*$xseq);
         $yin = pdl($r0*$yseq);
         plline($xin,$yin);
         $xout = pdl($r1*$xseq);
         $yout = pdl($r1*$yseq);
         plline($xout,$yout);
      }
   } elsif ($block_mode) {
      # loop over plot blocks
      for($i=0; $i<$num_plot_blocks; $i++) {
         $cur_blk = $index_good->index($i);
         # find boundaries of current block
         $x0 = $bnd_box->slice("0,0,$cur_blk")->sclr;
         $x1 = $bnd_box->slice("1,0,$cur_blk")->sclr;
         $y0 = $bnd_box->slice("0,1,$cur_blk")->sclr;
         $y1 = $bnd_box->slice("1,1,$cur_blk")->sclr;
         $x = pdl[$x0,$x0,$x1,$x1];
         $y = pdl[$y0,$y1,$y1,$y0];
         plline($x,$y);
      }
   }

   if ($vect_mode) {
      print "Plot vectors \n";
      $nrv = 30 if (! $nrv);
      $nav = 40 if (! $nav);
      $dr = $dist*($xmax-$xmin)/($nrv-1.);
      $dphi = 2.*pi/($nav);
      $velx_bin = zeroes $nrv,$nav;
      $vely_bin = zeroes $nrv,$nav;
      $rarr_bin = (sequence($nrv))->dummy(1,$nav);
      $rarr_bin = $dist*$xmin + $dr*$rarr_bin;
      $phiarr_bin = (sequence($nav))->dummy(0,$nrv);
      $phiarr_bin = -pi + $dphi*$phiarr_bin;
      # rebin the velocities
      rescale2d($dist*$temp_velx,$velx_bin);
#       $velx_bin = $temp_velx->map(t_identity, $velx_bin, {b=>'e'});
      rescale2d($dist*$temp_vely, $vely_bin);
#       $vely_bin = $temp_vely->map(t_identity, $vely_bin, {b=>'e'});
      if ($cart_mode) {
         $velxtmp = $velx_bin*cos($phiarr_bin) - $vely_bin*sin($phiarr_bin);
         $velytmp = $velx_bin*sin($phiarr_bin) + $vely_bin*cos($phiarr_bin);
         $cgrid2bin = plAlloc2dGrid ($rarr_bin * cos ($phiarr_bin),
            $rarr_bin * sin ($phiarr_bin));
         plvect ($velxtmp, $velytmp, 0.3, \&pltr2, $cgrid2bin);
      } else {
         $cgrid2bin = plAlloc2dGrid ($rarr_bin, $phiarr_bin);
         plvect ($velx_bin, $vely_bin , 0, \&pltr2, $cgrid2bin);
      }
   }

   if ($nbody) {
      $n_body_file = $filename;
      $n_body_file =~ s/hdf5/nbdy/;
      $n_body_file =~ s/cnt_//;
      read_nbody($n_body_file);
      if ($cart_mode) {
         plpoin($x_body,$y_body,1);
         plpoin(0,0,2);
      }
   }

   # $pl->shadeplot ($plot_var, $ns,
   #         BOX => ["$xrange[0]", "$xrange[1]", "$yrange[0]", "$yrange[1]"],
   #         PALETTE => 'RAINBOW', ZRANGE => ["$zmin","$zmax"],
   #         XBOX => 'BCNST');

   # $pl->xyplot (0,0, XBOX => 'BCTS', YBOX => 'BCTS');
   # plenv ("$xrange[0]", "$xrange[1]", "$yrange[0]", "$yrange[1]", 1, 0);
   # $pl->colorkey ($plot_var, 'v', VIEWPORT => [0.93, 0.96, 0.15, 0.85]);
   # $pl->colorkey ($shedge, 'h',
   #    VIEWPORT => [0.15, 0.85, 0.01, 0.03],
   #    XBOX => 'BC', YBOX => 'BCTS');

   # $pl->xyplot (0,0,
   #    viewport => [0.15, 0.85, 0.03, 0.05],
   #    xbox => 'bcts', ybox => 'bcs');

   # pladv(0);

   # draw color bar
#    if ($zmin != $zmax) {
   if (abs($zmin-$zmax) > abs($zmax)*1.e-8) {
      plvpor(0.1, 0.92, 0.03, 0.05);
      plwind($zmin, $zmax, 0.,1.);

      @box=($zmin, $zmax, 0., 1.);
      my $xinc = ($box[1] - $box[0])/$ncols;
      my $x0 = $box[0];

      for (my $i=0;$i<$ncols;$i++) {
         $x0 = $box[0] + ($i * $xinc);
         my $x1 = $x0 + $xinc;
         plcol1($i/$ncols);
         $x = pdl[$x0,$x0,$x1,$x1];
         $y = pdl[$box[2],$box[3],$box[3],$box[2]];
         plfill ($x,$y);
      }

      plcol0(1);
      plbox(0,0,0,0,"bctsm","bc");
   }
   # plfill(4,0.15, 0.85, 0.01, 0.03);

   # $pl->close;

   plend;
   plFree2dGrid ($cgrid2)

}

