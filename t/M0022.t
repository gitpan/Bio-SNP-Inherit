#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests=> 7;

BEGIN {
    use_ok('Bio::SNP::Inherit');
}

my $snp_obj_special = Bio::SNP::Inherit->new(
    manifest_filename => 't/M0022_F1_manifest.tab',
    data_filename     => 't/M0022_special_data.tab'
);

my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/M0022_F1_manifest.tab',
    data_filename     => 't/M0022_reduced_data.tab'
);

SKIP: {

    eval 'use File::Slurp qw{ slurp }';
    skip('File::Slurp needed for these tests',2) if $@;

    eval 'use Test::LongString';
    skip('Test::LongString needed for these tests',2) if $@;

    my $expected_special_abh = slurp('t/M0022_special_abh.tab');
    my $result_special_abh   = slurp('t/M0022_special_data.tab_abh.tab');
    is_string( $result_special_abh, $expected_special_abh,
        "ABH's corrected assigned" );

    my $expected_abh = slurp('t/M0022_reduced_abh.tab');
    my $result_abh   = slurp('t/M0022_reduced_data.tab_abh.tab');
    is_string( $result_abh, $expected_abh, "ABH's corrected assigned" );

}

#my $expected_summary = slurp ('t/M0022_summary.tab');
#my $result_summary = slurp('t/M0022_data.tab_summary.tab');
#is_string($result_summary, $expected_summary, 'summary correctly assigned');

_remove_temp_file('t/M0022_special_data.tab_summary.tab');
_remove_temp_file('t/M0022_special_data.tab_abh.tab');
_remove_temp_file('t/M0022_reduced_data.tab_summary.tab');
_remove_temp_file('t/M0022_reduced_data.tab_abh.tab');

sub _remove_temp_file {
    my $filename = shift;
    my $removed  = unlink $filename;
    ok( $removed, "removed tempoary file '$filename'" );
}
