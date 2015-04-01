<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html>
<HEAD>
	<META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html">
	<TITLE>HoloMask - Generate</TITLE>
</HEAD>
<body>


<?php
#HoloMask
#Computes a holographic mask for null testing mirrors
#Copyright Mauritz Andersson, 2004-2015
#Version: Dec 10, 2004
#Version: Dec 7, 2012
#Reseased as Open Source (MIT License): March 31, 2015

/*
The MIT License (MIT)

Copyright (c) [2015] [Mauritz Andersson]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/



#Get data from html post
$order=$_POST['order'];
#aperture radius
$diam=$_POST['diam'];
$r=$diam/2.;
#Radius of curvature
$R=$_POST['radi'];
#Conical const
$k=$_POST['k'];
#wavelength
$lamb=$_POST['lambda']*1e-6;
#Max slit size
$slit_size=$_POST['slit'];


#Computation units are mm
#convert to 1/72 inch (ps std)
$SCALE=72./25.4;
#center square half size
$center=0.3*$SCALE;

$ERRLIM=0.001; #mm
$STEP_SIZE=1.;  #mm
$phs=0.25; #Overall Phase shift of fringes



$ref_zon=0.75;
#reference zone
$l=2*pow($r,2)/$R/2*$ref_zon*(-$k);

$required_spacing=pow($r,2)/$R*(1-$ref_zon)*$r/$R*(-$k);
$a=$slit_size+abs($required_spacing);

$mean_period_on_mirror=$R/$a*$lamb;
$max_period_on_mirror=$R/$slit_size*$lamb;
$min_period_on_mirror=($R/($a+abs($required_spacing)))*$lamb;

$deviation=2*$k*pow(($r/sqrt(2)),4)  / (8*pow($R,3))/$lamb;

$interval=50; #mm, interval between support lines


#Output information
echo "<h1>Results</h1>\n"; 
echo "<table>";
echo "<tr><td>Diffraction order m=</td><td>$order </td></tr>\n";
echo "<tr><td>Mirror diameter (mm) </td><td>$diam </td></tr>\n";
echo "<tr><td>Mirror radius of curvature (mm) </td><td>$R</td></tr>\n"; 
echo "<tr><td>Conic  </td><td>$k<br/>\n"; 
echo sprintf("<tr><td>Wavelength (nm) </td><td>%f</td></tr>\n",$lamb*1E6); 

echo "<tr><td>Diffraction order separation (mm)</td><td>";
echo $a;
echo "</td></tr>\n";
echo "<tr><td>Max slit size (mm)</td><td>";
echo $slit_size;
echo "</td></tr>\n";
echo "<tr><td>Resolution on mirror,max period (mm) </td><td>";
echo $max_period_on_mirror;
echo "</td></tr>\n";
echo "<tr><td>min period (mm)</td><td>";
echo $min_period_on_mirror;
echo "</td></tr>\n";
echo "<tr><td>Deviation from sphere (lambda)</td><td>";
echo $deviation;
echo "</td></tr>\n";

if ($deviation<>0) {
echo "<tr><td>Wavelength and tolerance, PV=lambda/10 (nm)</td><td>";
echo $lamb*1E6;
echo " +- " ;
$dev=abs($lamb*(1.0/20.0/$deviation)*1E6);
echo $dev;
echo "</td></tr>\n";
}
else {
echo "<tr><td><b>Notice: You have a sphere and can null test without a mask.</b></td><td>";
echo "</td></tr>\n";
}

echo "</table>\n";

#Header for the EPS file
$header="%!PS-Adobe-2.0 EPSF-1.2
%%Creator:Mauritz Andersson, HoloMask
%%Title:HoloMask
%%CreationDate: 05/03/102 14:16
%%DocumentProcSets:Adobe_Illustrator_1.0 0 0
%%DocumentSuppliedProcSets:Adobe_Illustrator_1.0 0 0
%%DocumentFonts:Helvetica
%%BoundingBox: "
. sprintf("%d %d %d %d",0,0,round(2*$r*$SCALE)+1,round(2*$r*$SCALE)+1) .
"
%%TemplateBox:  "
. sprintf("%d %d %d %d",0,0,round(2*$r*$SCALE)+1,round(2*$r*$SCALE)+1) .
"
%%EndComments
%%BeginProcSet:Adobe_Illustrator_1.0 0 0
% Copyright (C) 1987 Adobe Systems Incorporated.
% All Rights Reserved.
% Adobe Illustrator is a trademark of Adobe Systems Incorporated
/Adobe_Illustrator_1.0 dup 100 dict def load begin
/Version 0 def
/Revision 0 def
% definition operators
/bdef {bind def} bind def
/ldef {load def} bdef
/xdef {exch def} bdef
% graphic state operators
/_K {3 index add neg dup 0 lt {pop 0} if 3 1 roll} bdef
/_k /setmybcolor where {/setcmybcolor get} 
	{{1 sub 4 1 roll _K _K _K setrgbcolor pop} bind} ifelse def
/g {/_b xdef /p {_b setgray} def} bdef
/G {/_B xdef /P {_B setgray} def} bdef
/k {/_b xdef /_y xdef /_m xdef /_c xdef /p {_c _m _y _b _k} def} bdef
/K {/_B xdef /_Y xdef /_M xdef /_C xdef /P {_C _M _Y _B _K} def} bdef
/d /setdash ldef
/_i currentflat def
/i {dup 0 equ {pop _i} if setflat} bdef
/j /setlinejoin ldef
/J /setlinecap ldef
/M /setmiterlimit ldef
/w /setlinewidth ldef
% path construction operators
/_R {.25 sub round .25 add} bdef
/_r {transform _R exch _R exch itransform} bdef
/c {_r curveto} bdef
/C /c ldef
/v {currentpoint 6 2 roll _r curveto} bdef
/V /v ldef
/y {_r 2 copy curveto} bdef
/Y /y ldef
/l {_r lineto} bdef
/L /l ldef
/m {_r moveto} bdef
% error operators
/_e [] def
/_E {_e length 0 ne {gsave 0 g 0 G 0 i 0 J 0 j 1 w 10 M [] 0 d
	/Courier 20 0 0 1 z [0.966 0.259 -0.259 0.966
	_e 0 get _e 2 get add 2 div _e 1 get _e 3 get add 2 div]
	e _f t T grestore} if} bdef
/_fill {{fill} stopped {/_e [pathbbox] def /_f 
	(ERROR: can't fill, increase flatness) def n _E} if} bdef
/_stroke {{stroke} stopped {/_e [pathbbox] def /_f 
	(ERROR: can't stroke, increase flatness) def n _E} if} bdef
% path painting operators
/n /newpath ldef
/N /n ldef
/F {p _fill} bdef
/f {closepath F} bdef
/S {P _stroke} bdef
/s {closepath S} bdef
/B {gsave F grestore S} bdef
/b {closepath B} bdef
% text block construction and painting operators
/_s /ashow /def
/_S {(?) exch {2 copy 0 exch put pop dup false charpath currentpoint 
	_g setmatrix _stroke _G setmatrix moveto 3 copy pop rmoveto} forall 
	pop pop pop n} bdef
/_A {_a moveto _t exch 0 exch} bdef
/_L {0 _l neg translate _G currentmatrix pop} bdef
/_w {dup stringwidth exch 3 -1 roll length 1 sub _t mul add exch} bdef
/_z [{0 0} {dup _w exch neg 2 div exch neg 2 div} 
	{dup _w exch neg exch neg}] bdef
/z {_z exch get /_a xdef /_t xdef /_l xdef exch findfont exch scalefont 
	setfont} bdef
/_g matrix def
/_G matrix def
/_D {_g currentmatrix pop gsave concat _G currentmatrix pop} bdef
/e {_D p /t {_A _s _L} def} bdef
/r {_D P /t {_A _S _L} def} bdef
/a {_D /t {dup p _A _s P _A _S _L} def} bdef
/o {_D /t {pop _L} def} bdef
/T {grestore} bdef
% group construction operators
/u {} bdef
/U {} bdef
% font construction operators
/Z {{findfont begin 
		currentdict dup length dict begin 
			1 index /FID ne {def} {pop pop} ifelse
	} forall 
	FontName exch def dup length 0 ne 
	{/Encoding Encoding 256 array copy def 0 exch {dup type /nametype eq
	{/Encoding 2 index 2 index put pop 1 add} {exch pop} ifelse} forall} if pop
	currentdict dup end end /FontName get exch definefont pop} bdef
	end
%%EndProcSet
%%EndProlog
%%BeginSetup
Adobe_Illustrator_1.0 begin
n
%%EndSetup
0 G
0 g
0 j
0 J
0 w	
0.000 0.000 0.000 1.000 k
";

#Writing file
$uniq=time();
$filename="data/holomask_$uniq.eps";
$handle = fopen($filename, "w");


fwrite($handle, $header);


#Start computing holomask



function z($x,$y)
{
#defines the wavefront
    global $a,$R,$k,$l,$lamb,$order;
    $z=(1/$order)*( $x*$a/$R + $l*(pow($x,2)+pow($y,2))/(2*pow($R,2)) + 2*$k*pow(pow($x,2)+pow($y,2),2) / (8*pow($R,3)) )/ $lamb ;
    return $z;
}


function z_x($x,$y)
{
#defines the x differentiation of the wavefront
    global $a,$R,$k,$l,$lamb,$order;
    $z_x=(1/$order)*( $a/$R + $l*(2*$x)/(2*pow($R,2)) + 2*$k*2*(pow($x,2)+pow($y,2))*2*$x  / (8*pow($R,3)) )/ $lamb; 
    return $z_x;
}


function find_x($x,$y,$phs)
{
    global $ERRLIM;
    #solve iteratively to find x pos of a fringe given y pos
    #using x as initial value
    $err=1e10;
    $count=0;
    while ( ($err>$ERRLIM) & ($count<60) )
    {
        $count=$count+1;
        $x_n=$x-(z($x,$y)-$phs)/z_x($x,$y);
        $err=abs($x_n-$x);
        $x=$x_n;
    }
    if ( abs(z($x,$y)-$phs)<(10*$ERRLIM))
    {
        return $x;
    } else
    {
#        echo "Holomask: Problem to find fringe.";
        return 1e20;
    }
}




function curve($x,$phs,$y_step)
{
    global $r;
    #create fringe segment at phs
    #x initial point (x need only be a guess)
    $i=0; 
    $upp=TRUE;
    #initial moveto point
    $y=$y_step/2.;
    $xx=find_x($x,$y,$phs);
    $x_pre=$xx;
    $y_pre=$y;

    $pts=array();
    array_push($pts,array($xx,$y) );
    while ($upp) {
        $i=$i+1;
        $yy=$y+$i*$y_step;
        $xx=find_x($x_pre,$yy,$phs);
        if ((pow($xx,2)+pow($yy,2))>pow($r,2)) {
            $upp=FALSE;
            #outside aperture, find intersection
            $kk=($xx-$x_pre)/($yy-$y_pre);
            $p=(2*$x_pre*$kk-2*$y_pre*pow($kk,2))/(1+pow($kk,2));
            $q=(pow($y_pre,2)*pow($kk,2)+pow($x_pre,2)-2*$x_pre*$y_pre*$kk-pow($r,2))/(1+pow($kk,2));
            $yy=-$p/2.+sqrt(max(pow($p/2.,2)-$q,0.));
            $xx=find_x($x_pre,$yy,$phs);
        }
        array_push($pts,array($xx,$yy) );
        $x_pre=$xx;
        $y_pre=$yy;
    }
    return $pts;
}


function path(&$points_l,&$points_r,$support_interval)
{
    #make path from two curves and remove elements
    unset($pth);
    if ( (count($points_l)>=$support_interval) & (count($points_r)>=$support_interval) ){
        $pth=array_splice($points_l,0,$support_interval);
        $tmp=array_splice($points_r,0,$support_interval);
        $tmp=array_reverse($tmp);
        $pth=array_merge($pth,$tmp);
    }
    else {
        $pth=array_splice($points_l,0,count($points_l));
        $tmp=array_splice($points_r,0,count($points_r));
        $tmp=array_reverse($tmp);
        $pth=array_merge($pth,$tmp);
    }
    return $pth;
}

function sign($s)
{
	return ($s)/abs($s+1e-16);
}

function fringe($i)
{
#Does the actual output
    global $r,$k,$phs,$SCALE,$STEP_SIZE,$interval,$handle;
    $more_path=TRUE;
    $crv_l=curve(-0.8*$r*sign($k),$i+$phs,$STEP_SIZE);
    $crv_r=curve(-0.8*$r*sign($k),$i+0.5+$phs,$STEP_SIZE);
    while ($more_path) {
        unset($pth);
        if (min(count($crv_l),count($crv_r))<(1.1*$interval)) {
            $pth=path($crv_l,$crv_r,1e10);
        } else {            
            $pth=path($crv_l,$crv_r,$interval);
        }
        if (count($pth)>2) {
                #up
                $ptm=array_pop($pth);
                $xxm=$ptm[0];
                $yym=$ptm[1];
                #moveto
                $str=sprintf("%01.2f %01.2f m\r\n",($xxm+$r)*$SCALE,($yym+$r)*$SCALE);
                fwrite($handle,$str);
                foreach ( $pth as $pt) {
                   $xx=$pt[0];
                   $yy=$pt[1];
                   #lineto
                   $str=sprintf("%01.2f %01.2f l\r\n",($xx+$r)*$SCALE,($yy+$r)*$SCALE);
                   fwrite($handle,$str);
                }
                #close path
                fwrite($handle,"s\r\n");
                #down
                #moveto
                $str=sprintf("%01.2f %01.2f m\r\n",($xxm+$r)*$SCALE,(-$yym+$r)*$SCALE);
                fwrite($handle,$str);
                foreach ( $pth as $pt) {
                   $xx=$pt[0];
                   $yy=$pt[1];
                   #lineto
                   $str=sprintf("%01.2f %01.2f l\r\n",($xx+$r)*$SCALE,(-$yy+$r)*$SCALE);
                   fwrite($handle,$str);
                }
                #close path
                fwrite($handle,"s\r\n");
        } else {
            $more_path=FALSE;
        }
    }
}


function fringes()
{
    global $r,$phs,$min_period_on_mirror;
    $count=0;
    for ($i=-round(2*$r/$min_period_on_mirror);$i<=round(2*$r/$min_period_on_mirror);$i++) {
        $x_l=min(find_x(-0.9*$r,0.,$i-0.1+$phs),find_x(-0.9*$r,0.,$i+0.6+$phs));
        $x_r=max(find_x(0.8*$r,0.,$i+0.6+$phs),find_x(0.8*$r,0.,$i-0.1+$phs));
        if ( ($x_l>=-$r) & ($x_r<=$r) ) {
            fringe($i);
            $count=$count+1;
	}
    }
        
    echo "<tr><td>Number of fringes </td><td>$count </td></tr>\n";
}


#Make the call
fringes();

$centerspot=sprintf("%f %f m\r\n %f %f l\r\n %f %f l\r\n %f %f l\r\n f\r\n ",(-$center+$r*$SCALE),(-$center+$r*$SCALE),($center+$r*$SCALE),(-$center+$r*$SCALE),($center+$r*$SCALE),($center+$r*$SCALE),(-$center+$r*$SCALE),($center+$r*$SCALE),(-$center+$r*$SCALE),(-$center+$r*$SCALE));
fwrite($handle,$centerspot);
$cornerspot=sprintf("%f %f m\r\n %f %f l\r\n %f %f l\r\n  f\r\n ",(0),(0),($center*10),(0),(0),($center*10),(0),(0));
fwrite($handle,$cornerspot);



$footer="
%%Trailer
_E end
";

fwrite($handle, $footer);
fclose($handle);

echo "<h2>Get file <a href=$filename>here</a></h2>\n";


###############

#######################


?>
Fileformat is PostScript (EPS) which can be read using e.g. <a href="http://www.acrobat.com">Acrobat</a>, <a href="http://www.cs.wisc.edu/~ghost/gsview/"> GSview</a> or Photoshop.


</body>

</html>