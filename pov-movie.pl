#!/usr/bin/perl
use Getopt::Std;
use Sys::Hostname;

open A,"$ARGV[0]/zinfo";
#srand 1729;
getopts("b:de:hmn:q:s:vx:z:");

# Print help message if -h option is specified
if (defined($opt_h)) {
	print "Usage: pov-movie.pl <switches> <snapshot-directory> [<header_file_name>]\n\n";
	print "Switches:\n";
	print "-b <frame>        (Start rendering from this frame)\n";
	print "-d                (Don't duplicate existing files)\n";
	print "-e <num>          (Only render every <num> frame)\n";
	print "-h                (Print help information)\n";
	print "-m                (Switch on movie creation)\n";
	print "-n <frame>        (Render up to <frame> frames)\n";
	print "-q <quality>      (Quality of rendering, 1=good, 3=extreme)\n";
	print "-s <frame>        (Render a single frame)\n";
	print "-v                (Verbose output)\n";
	print "-x <nr>           (Render number <nr> of realizations)\n";
	print "-z <prefix>       (The prefix of theRender number <nr>)\n";
	exit 0;
}

die "One or two arguments required" unless @ARGV==1 || @ARGV==2;

# Set output-related variables
$dr=$ARGV[0];
$hn=@ARGV==2?$ARGV[1]:"grPr";
$ff=$opt_l?"li":($opt_a?"an":"fr");
$as=defined $opt_c?$opt_c:-1;
$povm="pov_headers/$hn.pov";
$prd_len = defined $opt_o?$opt_o:0;
$fpref = defined $opt_z?$opt_z : "f";
$fsuff = "_nr";

# Set POV-Ray-related variables
$pov_width=1000;
$pov_height=1000;
$pov_flags=$opt_q<=1?($opt_q==1?"+R3 +A0.01 -J":"+R2 +A0.3 -J")
		    :($opt_q==2?"+R6 +A0.001 -J":"+R9 +A0.0001 -J");
if($opt_q==3) {$pov_width=1400;$pov_height=1400;}
$verb=$opt_v?"":">/dev/null 2>/dev/null";

# Loop over the available frames
$ntrans=0;
$every=$opt_e?$opt_e:1;
$h=1;
#starting frame gets over written by $opt_s
if(defined $opt_s) {
    $opt_b = $opt_s;
    $opt_n = 1;
}
$a=defined $opt_b?$opt_b:0;
$nr=0;
$ncam=0;
$in_fn=sprintf("%s/%s.%05d%s%d", $dr, $fpref, $a, $fsuff, $nr);
while(-e $in_fn) {

	# Assemble output filename and check for skip conditions
	$fn=sprintf "${ff}_%05d.png",$a;
    if ($nr>0) {
        $fn = sprintf "${ff}_%05d_${nr}.png", $a;
    }
	last if defined $opt_n && $a>=$opt_b+$opt_n;
    if ((defined $opt_d && -e "$dr/$fn") || (defined $opt_b && $a<$opt_b)) {

        if(defined $opt_x) {
            if($nr == $opt_x-1) {
                $nr=0;
                $a+=$every;
            } else {
                $nr+=1;
            }
        } else {
            $a += $every;
        }
        next;
    }

	$mrsq=0;

	# Find camera size and read in bacteria
	open C,$in_fn;
	$i=0;$cam=8*8;
	while(<C>) {
		($t,$x,$y,$l,$theta,$id_,$pid_,$on_,$sheath_,$cfac_)=split;
		next if($t =~ /#/);
        if (!defined $opt_a || (defined $opt_a && ($t%2)==$opt_c)) {
            $xx=$l*cos $theta;
            $yy=$l*sin $theta;
            $gx[$i]=$x-$xx;$gy[$i]=$y-$yy;
            $hx[$i]=$x+$xx;$hy[$i]=$y+$yy;
            $ty[$i]=$t;
            $o=$gx[$i]*$gx[$i]+$gy[$i]*$gy[$i];$cam=$o if $cam<$o;
            $o=$hx[$i]*$hx[$i]+$hy[$i]*$hy[$i];$cam=$o if $cam<$o;
            $id[$i]=$id_, $parent[$id_]=$pid_ if $opt_l;
            $i++;
        }
	}
	close C;
	$cam=sqrt $cam;
	$cam+=1;
	$cam*=2.5;
	$ncam=$cam if $ncam<$cam;

	# Assemble the POV file
	open A,">$dr/rtemp$h.pov" or die "Can't open temporary POV file\n";
	open B,$povm or die "Can't open master POV file\n";
	while(<B>) {
		s/CAM/$ncam/g;
		if(/^#include "bac\.pov"$/) {
			foreach (0..($i-1)) {
				print A <<EOF;
union {
	sphere{<$gx[$_],$gy[$_],0>,r}
	cylinder{<$gx[$_],$gy[$_],0>,<$hx[$_],$hy[$_],0>,r}
	sphere{<$hx[$_],$hy[$_],0>,r}
EOF
				if($opt_l) {
					if($ntrans==0) {
						print A "\ttexture{pigment{rgb 0.9} finish{f0}}\n}\n";
					} else {
						$id_=$id[$_];
						$id_=$parent[$id_] while $id_>=$ntrans;
						print A "\ttexture{pigment{rgb <$col[$id_]>} finish{f0}}\n}\n";
					}
				} elsif($opt_a) {
					print A "\ttexture{pigment{rgb <$col[$th[$_]]>} finish{f0}}\n}\n";
				} elsif($opt_f) {
					if($ty[$_]%2==0) {print A "\ttexture{pigment{rgb <$colf1[$cfac[$_]]>} finish{f0}}\n}\n";}
					if($ty[$_]%2==1) {print A "\ttexture{pigment{rgb <$colf2[$cfac[$_]]>} finish{f0}}\n}\n";}
				} elsif($opt_g) {
					if($ty[$_]%2==0) {print A "\ttexture{pigment{rgb <$colg1[$sheath[$_]]>} finish{f0}}\n}\n";}
					if($ty[$_]%2==1) {print A "\ttexture{pigment{rgb <$colg2[$sheath[$_]]>} finish{f0}}\n}\n";}
                }else {
					print A "\ttexture{t$ty[$_]}\n}\n";
				}
			}
			next;
		}
		print A;
	}
	close B;
	if($prd_len >0){
	print A <<EOF;
#declare ccr=0.3;
#declare ccsr=0.3;
union{
	cylinder {<-$prd_len,-$prd_len,0>,<$prd_len,-$prd_len,0>,ccr}
	cylinder {<-$prd_len,$prd_len,0>,<$prd_len,$prd_len,0>,ccr}
	cylinder {<-$prd_len,-$prd_len,0>,<-$prd_len,$prd_len,0>,ccr}
	cylinder {<$prd_len,-$prd_len,0>,<$prd_len,$prd_len,0>,ccr}
	sphere {<-$prd_len,-$prd_len,0>,ccsr}
	sphere {<-$prd_len,$prd_len,0>,ccsr}
	sphere {<$prd_len,-$prd_len,0>,ccsr}
	sphere {<$prd_len,$prd_len,0>,ccsr}
	texture{
		pigment{rgb <1,1,1> }
		finish{f0}
	}
}
EOF
	}
	close A;

    print "Rendering frame $a nr $nr...\n";
    exec "cd $dr; povray -D +O$fn +W$pov_width +H$pov_height $pov_flags rtemp$h.pov $verb" if (($pid[$h]=fork)==0);

    if(defined $opt_x) {
            if($nr == $opt_x-1) {
                $nr=0;
                $a+=$every;
            } else {
                $nr+=1;
            }
    } else {
        $a += $every;
    }
    # update the filename
    $in_fn=sprintf("%s/%s.%05d%s%d", $dr, $fpref, $a, $fsuff, $nr);
}
print "Done with rendering all frames.\n";

# Automatically create a movie if requested
if ($opt_m) {
	$o=$dr;
	$o=~s/\.out//;
    system "ffmpeg -r 20 -y -i $dr/${ff}_%05d.png -preset veryslow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $o.mov";
}

# Rainbow color scheme
sub rainbow2 {
	#print "@_[0]\n";
	$re=0.5+0.5*cos(@_[0]);
	$gr=0.5+0.5*cos(@_[0]-2.09439510239319549230842892219);
	$bl=0.5+0.5*cos(@_[0]-4.18879020478639098461685784436);
	return sprintf "%5.3f,%5.3f,%5.3f",rfunc($re),rfunc($gr),rfunc($bl);
}

sub rfunc {
	return @_[0]*@_[0]*(3-2*@_[0]);
}

sub gradient1{
	$re= 0.95+@_[0]*0.05;
	$gr= @_[0];
	$bl= @_[0];
	return sprintf "%6.4f,%6.4f,%6.4f",$re,$gr,$bl;
}

sub gradient2{
	$re= @_[0];
	$gr= @_[0];
	$bl= 0.95+0.05*@_[0];
	return sprintf "%6.4f,%6.4f,%6.4f",$re,$gr,$bl;
}
