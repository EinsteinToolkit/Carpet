#! /usr/bin/perl -s
#
# Blame Ian Hawke for not checking that Scott had already written a C
# program to do this add so writing a perl version instead
#
# Given an output file in CarpetIOASCII 1d format, strip to ygraph format.
# The arguments should be direction (x=0,y=1,z=2) and filename.
# Output is to a file. Only the base name should be given.
# The base will be appended with _<level>.xg
#
# Example:
#
# Carpet2ygraph.pl 0 alp.xl alp_x
#
# produces alp_x_<d>.xg where <d> is the refinement level.
#
# The headers of the Carpet files are just copied with the addition of
# the correct time for ygraph, and the different comment marker.
#
# This script is a much altered version of Tom Goodale's convergence
# testing scripts.
#

use FileHandle;

if(@ARGV != 3)
{
  print "Usage: $0 direction <Inputfile> <Outputfile> \n";
  exit;
}

open(CARPETFILE,   "<$ARGV[1]") || die "Unable to find file \"$ARGV[1]\"."; 

#
# Open the output file for the base grid. 
#

$file = $ARGV[2]."_0.xg";
my $fh = new FileHandle(">$file") || die "Unable to open file \"$file\"."; 
push(@outputfilelist,$fh);

#
# Find the correct column for the spatial coordinate; requires a magic number
#

$direction = $ARGV[0]+9; 

$flag = 0;
$timeset = 0;
$refinementlevel = 0;
$componentflag[0] = 0;
while (<CARPETFILE>)
{
    $line = $_;

    if(/^\#/) # The line is a header comment
    {
	if ($flag==1) # It's a new level and there is data to output
        {
	    $fh = $outputfilelist[$refinementlevel];
	    print $fh @outputdata;
	    @outputdata=("\n");
	    $flag = 0;
	}
	if ($line =~ /refinement level ([0-9])/) # Line gives ref. level
	{
	    $refinementlevel = $1;
	    $line =~ /component ([0-9+])/;
	    $componentflag[$refinementlevel] = $1; # Which component?

	    #
	    # If no file exists for this refinement level, 
	    # open and add filehandle to array
	    #

	    if ($refinementlevel > $#outputfilelist)
	    { 
		for ($i = $#outputfileflist+1; $i < $refinementlevel+1; $i++)
		{
		    $file = $ARGV[2]."_".$i.".xg";
		    my $fh = new FileHandle(">$file") || die "Unable to open file \"$file\"."; 
		    $outputfilelist[$i]=$fh;
		}
	    }
	}
	# Only output the headers if this is the zero component
	# FIXME: what happens if component 0 isn't output first?
	if (0 == $componentflag[$refinementlevel])
	{
	    push(@outputdata, ("\"",$line)); # Add ygraph comment marker
	}
	else
	{
	    $flag = 1;
	    @outputdata=("");
	}
    }
    else # The line contains real data
    {
	@data = split([/ \t]+/,$line);
	if ($flag== 0) # This is the first line of data
	{
	    $flag = 1;
	    $timeset = $data[8]; # Magic number gives the Cactus time
	    @outputdata = ("\n\"Time = ",$timeset,@outputdata);
	}
	push(@outputdata, $data[$direction], " ", $data[12]);
    }
}

#
# At end of file, so output final data set.
#

$fh = $outputfilelist[$refinementlevel];
print $fh @outputdata;

foreach $fh (@outputfilelist)
{
    close($fh);
}
