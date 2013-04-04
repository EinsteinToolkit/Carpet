#!/usr/bin/perl -w
use strict;

# Thorn Requirements can produce debug output describing which regions
# of which variables are valid and which are not at what time may be
# confusing. This is enabled via setting the parameter
# Requirements::output_changes = yes.
#
# This generates lines, prefixed with "LOC: " that contain information
# what is changed when and why. While this is sort of human-readable
# it's not ideal either. Because of that I created a script which
# takes that (just pipe the Cactus stdout into it) and spits out an
# html web page which displays the same information a bit nicer. An
# example can be seen here:
# 
#    https://www.cct.lsu.edu/~knarf/Requirements_example.html
# 
# (green is for a changed setting, bright color is for undefined, dark
# for defined)
# 
# Obviously this can be improved upon. Adding the name of the variable
# (instead of only the index (vi) would be an example. Thinking of
# something for the case of a gazillion GFs would also be needed.



sub col
{
  my ($changed, $to) = @_;
  if ($changed =~ /[a-z]/) {
    if ($to eq "0") {
      return "#CCCCCC";
    } else {
      return "#777777";
    }
  } else {
    if ($to eq "0") {
      return "#00CC00";
    } else {
      return "#007700";
    }
  }
}

my @lines = <>;

my $max_reflevel = 0;
foreach (@lines) {
  if ($_ =~ /LOC: vi:(-?\d+),it:(-?\d+), \[rl:(-?\d+),tl:(-?\d+),m:(-?\d+)\] \((in|IN):([01]),(bo|BO):([01]),(gh|GH):([01]),(bg|BG):([01])\)(.*)/) {
    if ($3 > $max_reflevel) { $max_reflevel = $3; };
  }
}

print "<html><body><table border=1>\n";
print "<tr><th>it</th><th>vi</th><th>rl</th><th>tl</th><th>m</th><th>&nbsp;</th><th>&nbsp;</th></tr>\n";
foreach (@lines) {
  if ($_ =~ /LOC: vi:(-?\d+),it:(-?\d+), \[rl:(-?\d+),tl:(-?\d+),m:(-?\d+)\] \((in|IN):([01]),(bo|BO):([01]),(gh|GH):([01]),(bg|BG):([01])\)(.*)/) {
    my $vi = $1;
    my $it = $2;
    my $rl = $3;
    my $tl = $4;
    my $m  = $5;
    my $IN = $6;
    my $in = $7;
    my $BO = $8;
    my $bo = $9;
    my $GH = $10;
    my $gh = $11;
    my $BG = $12;
    my $bg = $13;
    my $info = $14;
    print "<tr><td>$it</td><td>$vi</td>\n";
    print "    <td>$rl</td><td>$tl</td><td>$m</td>\n";
    print "<td>\n";
    my $empty_s = "<td height=10px width=74px>&nbsp;</td>";
    my $empty_l = "<td height=20px width=74px>&nbsp;</td>";
    my $prefix_s  = $empty_s x ($rl);
    my $prefix_l  = $empty_l x ($rl);
    print "<table cellpadding=0 cellspacing=1px><tr>";
    print $prefix_s;
    print "<td height=10px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=20px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=10px bgcolor=".&col($BG, $bg).">&nbsp;</td>";
    print "<td height=10px width=20px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "</tr><tr>";
    print $prefix_l;
    print "<td height=20px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=20px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=20px width=10px bgcolor=".&col($GH, $gh).">&nbsp;</td>";
    print "<td height=20px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=20px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "</tr><tr>";
    print $prefix_s;
    print "<td height=1px width=10px bgcolor=".&col($BG, $bg).">&nbsp;</td>";
    print "<td height=1px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=1px width=10px bgcolor=".&col($GH, $gh).">&nbsp;</td>";
    print "<td height=1px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=1px width=10px bgcolor=".&col($BG, $bg).">&nbsp;</td>";
    print "</tr><tr>";
    print $prefix_l;
    print "<td height=20px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=20px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=20px width=10px bgcolor=".&col($GH, $gh).">&nbsp;</td>";
    print "<td height=20px width=20px bgcolor=".&col($IN, $in).">&nbsp;</td>";
    print "<td height=20px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "</tr><tr>";
    print $prefix_s;
    print "<td height=10px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=20px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=10px bgcolor=".&col($BG, $bg).">&nbsp;</td>";
    print "<td height=10px width=20px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "<td height=10px width=10px bgcolor=".&col($BO, $bo).">&nbsp;</td>";
    print "</tr></table></td>\n";
    print "<td>$info</td>\n";
    print "</tr>\n";
  }
}
print "</table></body></html>\n";
