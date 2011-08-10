#use strict;
use IO::File;
use Bio::SeqIO;
use Carp;


my($prefix) = "bar_";
my @filehandles;


# (1) quit unless we have the correct number of command-line args
my($num_args) = $#ARGV + 1;
if ($num_args != 2) {
  print "\nUsage: name.pl first_name last_name\n";
  exit;
}

# (2) we got two command line args, so assume they are the
# first name and last name
print @ARGV;
print "\n";
my($seqfile, $barfile) = @ARGV;
print "seqfile=".$seqfile;
print "\n";


my $seq;
my $bar;
my $bc;
my $element;

load_barcode_file ("../barcodes");

create_output_files;

my $seqio = Bio::SeqIO->new(-file => "<$seqfile", '-format' => 'Fastq') 
            or die "could not open $seqfile: $!\n";;
my $bario = Bio::SeqIO->new(-file => "<$barfile", '-format' => 'Fastq') 
            or die "could not open $barfile: $!\n";
while ($seq = $seqio->next_seq and $bar = $bario->next_seq) {
    my $bc = $bar->seq;
    print "barcode = ".$bc.":\n";

    # find the matching barcode; return the index in the list or undef if not found.

    print STDERR "sequence $bc_bases: \n" if $debug;

    my $best_barcode_mismatches_count = $barcodes_length;
    my $best_barcode_ident = undef;

    #Try all barcodes, find the one with the lowest mismatch count
    foreach my $barcoderef (@barcodes) {
	my ($ident, $barcode) = @{$barcoderef};

	# Get DNA fragment (in the length of the barcodes)
	# The barcode will be tested only against this fragment
	# (no point in testing the barcode against the whole sequence)
	my $sequence_fragment;
	$sequence_fragment = substr $bar->seq, 0, $barcodes_length;

	my $mm = mismatch_count($bar->seq, $barcode) ; 

	# if this is a partial match, add the non-overlap as a mismatch
	# (partial barcodes are shorter than the length of the original barcodes)
	$mm += ($barcodes_length - length($barcode)); 

	if ( $mm < $best_barcode_mismatches_count ) {
	    $best_barcode_mismatches_count = $mm ;
	    $best_barcode_ident = $ident ;
	}
    }

    $best_barcode_ident = 'unmatched' 
	if ( (!defined $best_barcode_ident) || $best_barcode_mismatches_count>$allowed_mismatches) ;
    
    print STDERR "sequence $seq_bases matched barcode: $best_barcode_ident\n" if $debug;

    $counts{$best_barcode_ident}++;
    
    #get the file associated with the matched barcode.
    #(note: there's also a file associated with 'unmatched' barcode)
    my $outfile = $files{$best_barcode_ident};

    print "outfile $outfile\n";
    $outfile->write_seq($seq);
}



#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
#
# Code lifted from fastx_barcode_splitter.pl
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

#
# Read the barcode file
#
sub load_barcode_file ($) {
	my $filename = shift or croak "Missing barcode file name";

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
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length($barcode)<=$allowed_mismatches;

		$barcodes_length = length($barcode) unless defined $barcodes_length;
		die "Error: found barcodes in different lengths. this feature is not supported yet.\n" 
			unless $barcodes_length == length($barcode);

	 	push @barcodes, [$ident, $barcode];

		if ($allow_partial_overlap>0) {
			foreach my $i (1 .. $allow_partial_overlap) {
				substr $barcode, ($barcodes_at_bol)?0:-1, 1, '';
	 			push @barcodes, [$ident, $barcode];
			}
		}
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


# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	$barcodes{'unmatched'} = 1 ;

	foreach my $ident (keys %barcodes) {
		my $new_filename = $newfile_prefix . $ident . $newfile_suffix; 
		$filenames{$ident} = $new_filename;
		my $file = Bio::SeqIO->new(-file => ">$outfile", '-format' => 'Fastq') 
		           or die "could not open $outfile: $!\n";
		$files{$ident} = $file ;
	}
}
