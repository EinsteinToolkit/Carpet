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

while (<CARPETFILE>)
{
    $line = $_;
    if ($line =~ /^#/)
    {
	#Do nothing for headers!
    }
    elsif ($line =~ /([0-9+-ed.]+)*/)
    {
	@dataline = split(/ +/,$line);
	if ($dataline[8] - $time > 1.e-10)
	{
	    if ($new)
	    {
		print $fh "\n\n\"Time = ",$time,"\n";
		my @sortedcoords = sort numerically (keys %data);
		my $localcoord;
		foreach $localcoord (@sortedcoords)
		{
		    print $fh $localcoord." ".$data{$localcoord};
		}
	    }
	    $new++;
	    $time = $dataline[8];

	}
	my $coord = $dataline[9+$ARGV[0]];
	my $val = $dataline[12];
	$data{$coord} = $val;
    }
}
print $fh "\n\"Time = ",$time,"\n";
my @sortedcoords = sort numerically (keys %data);
my $localcoord;
foreach $localcoord (@sortedcoords)
{
    print $fh $localcoord." ".$data{$localcoord};
}

sub numerically {$a <=> $b;}
