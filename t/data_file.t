#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests=> 3;

BEGIN{ use_ok('Bio::SNP::Inherit'); }

eval{
my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/manifest_all.tab',
    data_filename => 't/data_empty_ids.tab',
);
};

like($@, qr/data .* tag/xmsi, 'fail on no [Data] tag in file'); 


_remove_temp_file('t/data_empty_ids.tab_summary.tab');

sub _remove_temp_file {
    my $filename = shift;
    my $removed  = unlink $filename;
    ok( $removed, "removed tempoary file '$filename'" );
}
