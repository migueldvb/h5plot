#!/usr/bin/perl -w
# Plot initial velocities in shearingbox calculation

use PDL;
use PDL::Graphics::PLplot;

use constant nx => 20;
use constant ny => 12;

$x = (-5. + 10.*sequence(nx)/(nx-1))->dummy(1,ny);
$y = (-5. + 10.*sequence(ny)/(ny-1))->dummy(0,nx);


$r3 = sqrt($x*$x + $y*$y);
$r3 = $r3*$r3*$r3;

$vely = -$x*1.5*($r3-1.)/($r3+0.75);
$velx = -$y*6./($r3+3.);

print $x,$y,$r3,$velx,$vely;

$cgrid2bin = plAlloc2dGrid ($x, $y);

plinit;
plenv(-5.,5.,-5.,5.,1,0);
plvect($velx,$vely,0,\&pltr2,$cgrid2bin);
plend;
