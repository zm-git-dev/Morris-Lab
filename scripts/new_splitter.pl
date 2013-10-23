#!/usr/bin/perl
use strict;
use warnings;

use IO::File;
use Carp;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw($GzipError);


my $barcodes = "barcodes";
my $barcodes_length;
my @barcodes;
my $mismatches = 0;
my $debug = 0;
my $help;

my %counts;
my $seqio;
my $bario;
my %filenames;		# output file names
my %files;		# output file handles

sub usage {
  print STDERR "Unknown option: @_\n" if ( @_ );
  print STDERR "usage: program [--barcodes|-b FILENAME] [--mismatches|-m NUMBER] [--help|-?] [--debug|-d] <sequence_fastq> <index_fastq>\n";
  print STDERR "sequence and index files may be compressed with gzip";
  exit;
}

#
# Read the barcode file
#
sub load_barcode_file ($) {
    my $filename = shift or croak "Missing barcode file name";
    
    print "using barcode file \"", $filename, "\"\n" if $debug;
    open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
    while (<BCFILE>) {
	next if m/^#/;
	chomp;
	my ($ident, $barcode) = split ;
	
	$barcode = uc($barcode);
	
	# Sanity checks on the barcodes
	die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
	die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
	    unless $barcode =~ m/^[AGCT]+$/;

	die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n" 
	    unless $ident =~ m/^\w+$/;

	die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
	    "mismatches ($mismatches). This makes no sense. Specify fewer  mismatches.\n" 
	    if length($barcode)<=$mismatches;

	$barcodes_length = length($barcode) unless defined $barcodes_length;
	die "Error: found barcodes in different lengths. this feature is not supported yet.\n" 
	    unless $barcodes_length == length($barcode);

	push @barcodes, [$ident, $barcode];
    }
    close BCFILE;
    
    if ($debug) {
	print STDERR "barcode\tsequence\n";
	foreach my $barcoderef (@barcodes) {
	    my ($ident, $seq) = @{$barcoderef};
	    print STDERR $ident,"\t", $seq ,"\n";
	}
    }
}		       

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
#
# Code lifted from fastx_barcode_splitter.pl
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }



#-- prints usage if no command line parameters are passed or there is an unknown
#   parameter or help option is passed
usage() if ( @ARGV < 1 or
          ! GetOptions('help|?' => \$help, 
		       'barcodes|b=s' => \$barcodes, 
		       'mismatches|m=i' => \$mismatches, 
		       'debug|d' => \$debug)
          or defined $help );


# (1) quit unless we have the correct number of command-line args
my($num_args) = $#ARGV + 1;
if ($num_args != 2) {
  usage()
}

# (2) we got two command line args, so assume they are the
# first name and last name

my($seqfile, $barfile) = @ARGV;

if ($seqfile =~ /.*\.gz/) {
    print "uncompressing ", $seqfile, "\n" if $debug;
    $seqio = IO::Uncompress::Gunzip->new( $seqfile )
	or die "Cannot uncompress $seqfile: $GunzipError\n";
} else {
    print "opening ", $seqfile, "\n" if $debug;
    open( $seqio, "<", $seqfile ) || die "Can't open $seqfile: $!";
}

if ($barfile =~ /.*\.gz/) {
    print "uncompressing ", $barfile, "\n" if $debug;
    $bario = IO::Uncompress::Gunzip->new( $barfile )
	or die "Cannot uncompress $barfile: $GunzipError\n";
} else {
    print "opening ", $barfile, "\n" if $debug;
    open( $bario, "<", $barfile ) || die "Can't open $barfile: $!";
}

print $mismatches," allowed mismatches\n" if $debug;

load_barcode_file ($barcodes);

create_output_files();

while (my $bar = <$bario>) {
    my @bar;
    my @seq;

    $bar[0] = $bar;
    $bar[1] = <$bario>;
    $bar[2] = <$bario>;
    $bar[3] = <$bario>;

    my $best_barcode_mismatches_count = $barcodes_length;
    my $best_barcode_ident = undef;

    #Try all barcodes, find the one with the lowest mismatch count
    foreach my $barcoderef (@barcodes) {
	my ($ident, $barcode) = @{$barcoderef};

	# Get DNA fragment (in the length of the barcodes)
	# The barcode will be tested only against this fragment
	# (no point in testing the barcode against the whole sequence)
	my $sequence_fragment;
	$sequence_fragment = substr $bar[1], 0, $barcodes_length;

	my $mm = mismatch_count($sequence_fragment, $barcode) ; 

	# print "$mm comparing $sequence_fragment to $barcode\n";

	# if this is a partial match, add the non-overlap as a mismatch
	# (partial barcodes are shorter than the length of the original barcodes)
	$mm += ($barcodes_length - length($barcode)); 

	if ( $mm < $best_barcode_mismatches_count ) {
	    $best_barcode_mismatches_count = $mm ;
	    $best_barcode_ident = $ident ;
	}
    }

    $best_barcode_ident = 'unmatched' 
	if ( (!defined $best_barcode_ident) || $best_barcode_mismatches_count>$mismatches) ;
    
    $counts{$best_barcode_ident}++;
    
    #get the file associated with the matched barcode.
    #(note: there's also a file associated with 'unmatched' barcode)
    my $out = $files{$best_barcode_ident};

    # now read forward in the sequence file to find the corresponding sequence.

    my $gapcount = 0;
    my $spacer = "";
    my $barid = substr $bar, 0, -3;
    while (my $seq = <$seqio>) {
	$seq[0] = $seq;
	$seq[1] = <$seqio>;
	$seq[2] = <$seqio>;
	$seq[3] = <$seqio>;

	my $seqid = substr $seq, 0, -3;
	if ($seqid eq $barid) {
	    print $out $seq[0];
	    print $out $seq[1];
	    print $out $seq[2];
	    print $out $seq[3];
	    last;
	}
	++$gapcount;
	if (($gapcount % 100) == 0) {
	    print STDERR $spacer,$gapcount;
	    $spacer = "\r";
	}
    }
    if ($spacer ne "") {
	print STDERR "\n";
    }
    
}

close_output_files();


# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
    my $newfile_prefix = "";
    my $newfile_suffix = ".gz"; 

    #generate a uniq list of barcode identifiers;
    my %barcodes = map { $_->[0] => 1 } @barcodes; 
    $barcodes{'unmatched'} = 1 ;

    foreach my $ident (keys %barcodes) {
	my $new_filename = $newfile_prefix . $ident . $newfile_suffix; 
	$filenames{$ident} = $new_filename;
	print "creating output file for $filenames{$ident}\n" if $debug;
	my $file = new IO::Compress::Gzip $filenames{$ident}
	    or die "gzip failed to compress $filenames{$ident}: $GzipError\n";
	# open my $file, ">$filenames{$ident}" 
	#     or die "could not open $outfile: $!\n";
	$files{$ident} = $file ;
    }
}

# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub close_output_files {
    #generate a uniq list of barcode identifiers;
    my %barcodes = map { $_->[0] => 1 } @barcodes; 
    $barcodes{'unmatched'} = 1 ;

    foreach my $ident (keys %barcodes) {
	close($files{$ident});
    }
}


