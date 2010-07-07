#!/usr/bin/perl
use strict;
use warnings;
use Bio::SNP::Inherit;
if ( !@ARGV ) {
    print 'usage:' . "\n";
    print ' analyze.pl <manifest file name> <data file name>' . "\n";
}
else {
    my $snp_obj = Bio::SNP::Inherit->new(
        manifest_filename => shift,
        data_filename     => shift,
    );
}

=pod

=cut
