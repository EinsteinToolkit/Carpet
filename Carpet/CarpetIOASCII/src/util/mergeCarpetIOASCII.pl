#! /usr/bin/perl -w
#/*@@
#  @file      mergeCarpetIOASCII.pl
#  @date      Tue 15 August 2006
#  @author    Thomas Radke
#  @desc
#             Perl script to merge CarpetIOASCII output files,
#             eliminating duplicate datasets.
#
#             Source code comments also by Thomas Radke :-)
#  @enddesc
#@@*/


if (@ARGV) {
  print "\n" .
        "  This script can be used to merge CarpetIOASCII output written\n" .
        "  before and after recovery.\n" .
        "  It reads from STDIN the contents of one or more files\n" .
        "  in CarpetIOASCII format and writes them to STDOUT again,\n" .
        "  eliminating duplicate datasets.\n\n" .
        "  Example: cat alp.x.asc | $0 > alp.x.asc.merged\n\n";
  exit;
}

# Skalar-Variable zur Speicherung der aktuellen Zeilennummer
my $line = 0;

# Rauten-Feld zur Speicherung der bereits ausgegebenen Datensaetze
my %datasets = ();

# lese die naechste Zeile von der Standard-Eingabe, solange das Dateiende
# nicht erreicht ist
while (<STDIN>) {

  # aktualisiere die Zeilennummer
  $line++;

  # vergleiche die aktuelle Zeile mit der Markierungszeile
  # fuer einen neuen Datensatz
  if (/^# iteration (\d+)$/) {

    # ermittle die Iterationsnummer aus der Markierungszeile
    my $iteration = $1;

    # lies die naechste Zeile ein und pruefe, ob sie dem CarpetIOASCII-Format
    # entspricht
    $_ = <STDIN>;
    $line++;
    die "Format error in line $line: expected '# refinement level ...'\n"
      unless (/^# refinement level (\d+)   multigrid level (\d+)   map (\d+)   component (\d+)   time level (\d+)$/);

    # erzeuge einen eindeutigen Schluessel fuer den aktuellen Datensatz
    my $key = "$iteration $1 $2 $3 $4 $5";

    # pruefe im Rauten-Feld, ob dieser Datensatz bereits ausgegeben wurde
    while (defined $datasets{$key}) {

      # lies alle folgenden Zeilen von der Standard-Eingabe
      # bis zur Markierungszeile des naechsten Datensatzes
      while (<STDIN>) {
        $line++;
        next unless (/^# iteration (\d+)$/);

        # erzeuge einen eindeutigen Schluessel fuer diesen Datensatz
        $iteration = $1;
        $_ = <STDIN>;
        $line++;
        die "Format error in line $line: expected '# refinement level ...'\n"
          unless (/^# refinement level (\d+)   multigrid level (\d+)   map (\d+)   component (\d+)   time level (\d+)$/);
        $key = "$iteration $1 $2 $3 $4 $5";
        last;
      }
    }

    # gib die Markierungszeile des aktuellen Datensatzes aus
    print "# iteration $iteration\n";

    # vermerke im Rauten-Feld, dass dieser Datensatz schon ausgegeben wurde
    $datasets{$key} = 1;
  }

  # gib die aktuelle Zeile auf der Standard-Ausgabe aus
  print;
}
