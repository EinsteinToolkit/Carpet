#! /usr/bin/perl -sw
#
# Given an output file in CarpetIOScalar format or
# CarpetIOASCII scalar or 1d format, strip to ygraph format.

# Output is to a file named after the input file, with its extension ".asc"
# substituted by ".xg".
#
# Example:
#
#   Carpet2ygraphCat.pl alp.x.asc
#
# produces alp.x.xg
#
#   Carpet2ygraph.pl alp.minimum.asc
#
# produces alp.minumum.xg
#
# Some headers of the Carpet files are just copied with the addition of
# the correct time for ygraph, and the different comment marker.
#
# This script is an altered version of Ian Hawke's script.
#


use FileHandle;

if(@ARGV != 1)
{
  print "Usage: $0 <Inputfile>\n";
  exit;
}

my $direction;
my $filename = $ARGV[0];
open(CARPETFILE, "< $filename") || die "Unable to find input file '$filename'\n.";

if ($filename =~ /\.(x|y|z|d)\.asc$/) {
  $direction = ord($1) - ord('x') + 9;
  print "extracting CarpetIOASCII 1D data in direction $1...\n";
} elsif ($filename =~ /\.\.asc$/) {
  $direction = 8;
  print "extracting CarpetIOASCII 0D data...\n";
} elsif ($filename =~ /\.asc$/) {
  $direction = 0;
  print "extracting CarpetIOScalar data...\n";
} else {
  die "Do not recognize input file format for '$filename'.\n";
}

$filename =~ s/asc$/xg/;
my $fh = new FileHandle("> $filename") || die "Unable to open output file '$filename'.\n";

# deal with scalar output
if ($direction == 0) {
  while (<CARPETFILE>) {
    # The line is not a header comment
    if(/^[0-9]/) {
      my @data = split(/\s+/);
      print $fh $data[1]." ".$data[2]."\n";
    }
  }
  close ($fh);
  exit;
}


my %data;
my $time = -1;
my $new = 0;
my $currentit = -1;
my $lastit = -1;

my @datatoprint;
my $nsets = 1;
my $maxlength = 0;
my @lengths;

$lengths[0]=0;
while (<CARPETFILE>)
{
  chomp;
  next if (/^$/);

  if (/iteration/) {
    @itline = split(/ +/);
    $currentit = $itline[2];
  }
  #Do nothing for headers!
  next if (/^#/);

  @dataline = split(/[ \t]+/);
  if ($currentit != $lastit) {
    if ($new) {
      # do not print "Time..." for zero-D data
      push(@datatoprint,"\n\n\"Time = ".$time."\n") if ($direction !~ 8);

      my @sortedcoords = sort numerically (keys %data);
      foreach my $localcoord (@sortedcoords) {
        push(@datatoprint, $localcoord." ".$data{$localcoord}."\n");
      }
      $maxlength = $maxlength > (scalar @sortedcoords) ? $maxlength : (scalar @sortedcoords);
      $lengths[$nsets-1]=(scalar @sortedcoords);
      $nsets++;
      $lengths[$nsets-1]=0;
      %data=();
    }
    $new++;
    $time = $dataline[8];
    $lastit = $currentit;
  }
  my $coord = $dataline[$direction];
  my $val = $dataline[12];
  $data{$coord} = $val;
}

# do not print "Time..." for zero-D data
push(@datatoprint,"\n\n\"Time = ".$time."\n") if ($direction !~ 8);

my @sortedcoords = sort numerically (keys %data);
foreach my $localcoord (@sortedcoords) {
  push(@datatoprint, $localcoord." ".$data{$localcoord}."\n");
}
$maxlength = $maxlength > (scalar @sortedcoords) ? $maxlength : (scalar @sortedcoords);
$lengths[$nsets-1]=(scalar @sortedcoords);
$nsets++;
$lengths[$nsets-1]=0;

my $oldline="";
$nouts=0;
my $set=0;
foreach $line (@datatoprint) {
  if ($line =~ "Time") {
    if ($oldline) {
      for (my $i=$lengths[$set-1]; $i<$maxlength;$i++) {
        $nouts++;
        print $fh $oldline;
      }
    }
    $set++;
    print $fh $line;
  }
  else {
    $nouts++;
    print $fh $line;
    $oldline=$line
  }
}

# TR: I don't see why this is needed, it seems to print the last data item twice
# in CarpetIOBasic 0D output
#for (my $i=$lengths[$set-1]; $i<$maxlength;$i++) {
#  $nouts++;
#  print $fh $oldline;
#}

sub numerically {$a <=> $b;}
