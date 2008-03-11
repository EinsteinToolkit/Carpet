#! /usr/bin/perl -w

use strict;

sub draw ($$$$);

my $levels = 3;
my $xpos = 0;
my $ypos = 0;
my $size = 100 * 3**$levels;

my $linewidth = 10;
my $linecolour = 32;

# 1=green, 2=blue, 4=red, 5=magenta
my @colours = (4, 2, 1, 5);

my $x1 = $xpos;
my $y1 = $ypos;
my $x2 = $xpos + $size;
my $y2 = $ypos + $size;

print <<EOF;
\#FIG 3.2  Produced for xfig version 3.2.5
Landscape
Center
Inches
Letter  
100.00
Single
-2
1200 2
0 32 \#c0c0c0
2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5
	 $x1 $y1 $x2 $y1 $x2 $y2 $x1 $y2 $x1 $y1
EOF

draw 0, $xpos, $ypos, $size;
exit 0;



sub draw ($$$$)
{
    my ($level, $xpos, $ypos, $size) = @_;
    
    # coloured square
    {
        my $x1 = $xpos + $size/3;
        my $y1 = $ypos + $size/3;
        my $x2 = $xpos + 2*$size/3;
        my $y2 = $ypos + 2*$size/3;
        my $c = $colours[$level];
        
        print <<EOF;
2 2 0 1 $c $c 50 -1 20 0.000 0 0 -1 0 0 5
	 $x1 $y1 $x2 $y1 $x2 $y2 $x1 $y2 $x1 $y1
EOF
    }
    
    # horizontal grey line
    {
        my $x1 = $xpos;
        my $y1 = $ypos + $size/2 - $linewidth;
        my $x2 = $xpos + $size;
        my $y2 = $ypos + $size/2 + $linewidth;
        my $linecolour = 32;
        
        print <<EOF;
2 2 0 1 $linecolour $linecolour 100 -1 20 0.000 0 0 -1 0 0 5
	 $x1 $y1 $x2 $y1 $x2 $y2 $x1 $y2 $x1 $y1
EOF
    }
    
    # vertical grey line
    {
        my $x1 = $xpos + $size/2 - $linewidth;
        my $y1 = $ypos;
        my $x2 = $xpos + $size/2 + $linewidth;
        my $y2 = $ypos + $size;
        my $linecolour = 32;
        
        print <<EOF;
2 2 0 1 $linecolour $linecolour 100 -1 20 0.000 0 0 -1 0 0 5
	 $x1 $y1 $x2 $y1 $x2 $y2 $x1 $y2 $x1 $y1
EOF
    }
    
    # recur
    if ($level+1 < $levels) {
        for my $i (0, 1, 2) {
            for my $j (0, 1, 2) {
                if (! ($i==1 && $j==1)) {
                    draw $level+1, $xpos+$i*$size/3, $ypos+$j*$size/3, $size/3;
                }
            }
        }
    }
}
