use strict;
use warnings;
use Test::More tests=> 11;
BEGIN{
    use_ok('Bio::SNP::Inherit');
}

#Test that duplicate sample codes cause crash
eval{
        my $dead_snp_obj = Bio::SNP::Inherit->new(
            manifest_filename => 't/sample_bad_manifest.tab',
            data_filename => 't/sample_data.tab',
        );
};
like(
    $@,
    qr/duplicate sample id in manifest file/,
    'fail on duplicate sample id in manifest file'
);

#Test that having a missing id causes program to die
eval{
        my $dead_snp_obj = Bio::SNP::Inherit->new(
            manifest_filename => 't/missing_id_manifest.tab',
            data_filename => 't/sample_data.tab',
        );
};
like(
     $@,
     qr/Empty id/,
     'fail on missing id in manifest file'
);

my $dummy_sample_snp_obj = Bio::SNP::Inherit->new(
            manifest_filename => 't/dummy_manifest.tab',
            data_filename => 't/dummy_data.tab',
);

my %dummy_sample_for = (
    'WG0096796-DNAF01' => {
        name      => 'B97',
        group     => '25DL',
        replicates => ['WG0096796-DNAF11'],
    },
    'WG0096796-DNAA05' => {
        name      => 'B73xB97',
        group     => 'NAM F1',
        parentA    => 'WG0096795-DNAA01',
        parentB    => 'WG0096796-DNAF01',
    },
    'WG0096795-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicates =>['WG0096796-DNAA01','WG0096797-DNAA01'],
    },
    'WG0096796-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicate_of => 'WG0096795-DNAA01',
    },
    'WG0096797-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicate_of => 'WG0096795-DNAA01',
    },
    'WG0096796-DNAF11' => {
        name      => 'B97',
        group     => 'Holland pop',
        replicate_of => 'WG0096796-DNAF01',
    },
    'DUMMY-NUMBER-WEL' => {
        name      => 'B73xB97',
        group     => 'NAM F1',
        parentA    => 'WG0096795-DNAA01',
    },
    'DUMMY-NUMBE-WELL' => {
        name      => 'B73xB97',
        group     => 'NAM F1',
        parentB    => 'WG0096796-DNAF01',
    },
);

is_deeply(
            $dummy_sample_snp_obj->_sample_for(),
            \%dummy_sample_for,
            'correctly parses/stores info from dummy manifest file'
         );

SKIP: {

    eval 'use File::Slurp qw{ slurp }';
    skip('File::Slurp needed for these tests',2) if $@;

    eval 'use Test::LongString';
    skip('Test::LongString needed for these tests',2) if $@;

    my $expected_dummy_abh = slurp('t/dummy_abh.tab');
    my $result_dummy_abh   = slurp('t/dummy_data.tab_abh.tab');
    is_string( $result_dummy_abh, $expected_dummy_abh, 'dummy abh correct' );

    my $expected_dummy_summary = slurp('t/dummy_summary.tab');
    my $result_dummy_summary   = slurp('t/dummy_data.tab_summary.tab');
    is_string( $result_dummy_summary, $expected_dummy_summary,
        'dummy summary correct' );

}

#Test that information from the manifest file is correctly parsed and stored
#   in the object.
my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/sample_manifest.tab',
    data_filename     => 't/sample_data.tab',
);

my %expected_sample_for = (
    'WG0096796-DNAF01' => {
        name      => 'B97',
        group     => '25DL',
        replicates => ['WG0096796-DNAF11'],
    },
    'WG0096796-DNAA05' => {
        name      => 'B73xB97',
        group     => 'NAM F1',
        parentA    => 'WG0096795-DNAA01',
        parentB    => 'WG0096796-DNAF01',
    },
    'WG0096795-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicates =>['WG0096796-DNAA01','WG0096797-DNAA01'],
    },
    'WG0096796-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicate_of => 'WG0096795-DNAA01',
    },
    'WG0096797-DNAA01' => {
        name      => 'B73',
        group     => 'Control',
        replicate_of => 'WG0096795-DNAA01',
    },
    'WG0096796-DNAF11' => {
        name      => 'B97',
        group     => 'Holland pop',
        replicate_of => 'WG0096796-DNAF01',
    },
);

is_deeply(
            $snp_obj->_sample_for(),
            \%expected_sample_for,
            'correctly parses/stores info from sample manifest file'
         );

remove_temp_file('t/dummy_data.tab_abh.tab'     );
remove_temp_file('t/dummy_data.tab_summary.tab' );
remove_temp_file('t/sample_data.tab_abh.tab'    );
remove_temp_file('t/sample_data.tab_summary.tab');

sub remove_temp_file{
    my $filename = shift;
    my $result = unlink $filename;
    ok($result, "removed temp file: '$filename'");
}
