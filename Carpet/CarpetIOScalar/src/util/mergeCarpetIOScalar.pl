#! /usr/bin/perl -w
#/*@@
#  @file      mergeCarpetIOScalar.pl
#  @date      Tue 15 August 2006
#  @author    Thomas Radke
#  @desc
#             Perl script to merge CarpetIOScalar output files,
#             eliminating duplicate timesteps.
#
#             Source code comments also by Thomas Radke :-)
#  @enddesc
#@@*/


if (@ARGV) {
  print "\n" .
        "  This script can be used to merge CarpetIOScalar output written\n" .
        "  before and after recovery.\n" .
        "  It reads from STDIN the contents of one or more files\n" .
        "  in CarpetIOScalar format and writes them to STDOUT again,\n" .
        "  eliminating duplicate timesteps.\n\n" .
        "  Example: cat alp.norm1.asc | $0 > alp.norm1.asc.merged\n\n";
  exit;
}

# Skalar-Variable zur Speicherung der aktuellen Zeilennummer
my $line = 0;

# Rauten-Feld zur Speicherung der bereits ausgegebenen Datensaetze
my %timesteps = ();

# lese die naechste Zeile von der Standard-Eingabe, solange das Dateiende
# nicht erreicht ist
while (<STDIN>) {

  # aktualisiere die Zeilennummer
  $line++;

  # vergleiche die aktuelle Zeile mit dem Format fuer eine Datenzeile
  if (/^(\d+) /) {

    # ermittle die Iterationsnummer aus der Datenzeile
    my $iteration = $1;

    # pruefe im Rauten-Feld, ob dieser Datensatz bereits ausgegeben wurde
    while (defined $timesteps{$iteration}) {

      # lies alle folgenden Zeilen von der Standard-Eingabe
      # bis zur Datenzeile des naechsten Datensatzes
      while (<STDIN>) {
        $line++;
        next unless (/^(\d+) /);

        $iteration = $1;
        last;
      }
    }

    # vermerke im Rauten-Feld, dass dieser Datensatz ausgegeben wurde
    $timesteps{$iteration} = 1;
  }

  # gib die aktuelle Zeile auf der Standard-Ausgabe aus
  print;
}
