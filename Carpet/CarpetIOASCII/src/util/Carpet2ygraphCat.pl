#! /usr/bin/perl -s

use FileHandle;

if(@ARGV != 3)
{
  print "Usage: $0 direction <Inputfile> <Outputfile> \n";
  exit;
}

open(CARPETFILE,   "<$ARGV[1]") || die "Unable to find file \"$ARGV[1]\".";

$file = $ARGV[2].".xg";
my $fh = new FileHandle(">$file") || die "Unable to open file \"$file\".";

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
    $line = $_;
    if ($line =~ "iteration") {
      @itline = split(/ +/,$line);
      $currentit = @itline[2];
    }
    elsif ($line =~ /^#/)
    {
	#Do nothing for headers!
    }
    elsif ($line =~ /([0-9+-ed.]+)*/)
    {
	@dataline = split(/ +/,$line);
	if ($currentit != $lastit)
	{
	    if ($new)
	    {
		push(@datatoprint,"\n\n\"Time = ".$time."\n");
		my @sortedcoords = sort numerically (keys %data);
		my $localcoord;
		foreach $localcoord (@sortedcoords)
		{
		    push(@datatoprint, $localcoord." ".$data{$localcoord});
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
	my $coord = $dataline[9+$ARGV[0]];
	my $val = $dataline[12];
	$data{$coord} = $val;
    }
}
push(@datatoprint, "\n\"Time = ",$time,"\n");
my @sortedcoords = sort numerically (keys %data);
my $localcoord;
foreach $localcoord (@sortedcoords)
{
    push(@datatoprint, $localcoord." ".$data{$localcoord});
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
for (my $i=$lengths[$set-1]; $i<$maxlength;$i++) {
  $nouts++;
  print $fh $oldline;
}

sub numerically {$a <=> $b;}
