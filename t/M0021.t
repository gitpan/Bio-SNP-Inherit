#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests=> 5;

BEGIN {
    use_ok('Bio::SNP::Inherit');
}

my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/manifest_M0021.tab',
    data_filename     => 't/M0021_data.tab'
);

my @result_codes = $snp_obj->_ids_from_data_file();

SKIP: {

    eval 'use File::Slurp qw{ slurp }';
    skip('File::Slurp needed for these tests',2) if $@;

    eval 'use Test::LongString';
    skip('Test::LongString needed for these tests',2) if $@;

    my $expected_abh = slurp('t/M0021_abh.tab');
    my $result_abh   = slurp('t/M0021_data.tab_abh.tab');
    is_string( $result_abh, $expected_abh, "ABH's corrected assigned" );

    my $expected_summary = slurp('t/M0021_summary.tab');
    my $result_summary   = slurp('t/M0021_data.tab_summary.tab');
    is_string( $result_summary, $expected_summary,
        'summary correctly assigned' );

}

_remove_temp_file('t/M0021_data.tab_summary.tab');
_remove_temp_file('t/M0021_data.tab_abh.tab');

sub _remove_temp_file {
    my $filename = shift;
    my $removed  = unlink $filename;
    ok( $removed, "removed tempoary file '$filename'" );
}
