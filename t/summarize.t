use strict;
use warnings;
use Test::More tests=> 7;

BEGIN {
    use_ok('Bio::SNP::Inherit');
}

#Note that object creation causes processing of the files and the output of 
# a summary file and a file describing the samples' allelic inheritance.
my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/sample_manifest.tab',
    data_filename     => 't/sample_data.tab'
);




#INTERNAL METHODS (i.e. not part of the public interface, so they may change in the future.)

my @expected_ids = qw{
    WG0096795-DNAA01
    WG0096796-DNAA01
    WG0096797-DNAA01
    WG0096796-DNAA05
    WG0096796-DNAF01
    WG0096796-DNAF11};

my @result_ids = $snp_obj->_ids_from_data_file();

is_deeply( \@result_ids, \@expected_ids,
    'sample ids found in data file' );


my @columns = $snp_obj->_get_columns_for('WG0096795-DNAA01');
is_deeply( 
            \@columns, 
            [0,1,2], 
            'get columns for sample id '
            . '(i.e. original column of data plus columns of "replicates")'
         );


my $expected_summary_header = 
                "\t" 
                . join("\t", qw{ 
                AA
                AC
                AG
                AT
                CC
                CG
                CT
                GG
                GT
                TT
                --
                inconsistent
                parentA_polymorphic
                parentB_polymorphic
                parentA_unknown
                parentB_unknown
                sum_of_odds_and_ends
                F1_missing
                }) . "\n";
is( $snp_obj->_summary_header(), $expected_summary_header, 'summary header');


SKIP: {

    eval 'use File::Slurp qw{ slurp }';
    skip('File::Slurp needed for these tests',1) if $@;

    eval 'use Test::LongString';
    skip('Test::LongString needed for these tests',1) if $@;

#Test summary output file
my $expected_summary = slurp('t/sample_info.tab');
my $result_summary = slurp('t/sample_data.tab_summary.tab');
is_string($result_summary, $expected_summary, 'summary file created correctly');

}

_remove_temp_file('t/sample_data.tab_summary.tab');
_remove_temp_file('t/sample_data.tab_abh.tab');

sub _remove_temp_file {
    my $filename = shift;
    my $removed  = unlink $filename;
    ok( $removed, "removed tempoary file '$filename'" );
}

##Test allelic inheritance file
##TODO: Update abh file for new inheritance codes
##
#TODO: {
#    local $TODO = 'Update abh file for new inheritance codes' . "\n"
#        . "I'm not worrying with this at the moment, since M0022.t covers abh really well";
#
#SKIP: {
#
#        eval 'use File::Slurp qw{ slurp }';
#        skip('File::Slurp needed for these tests') if $@;
#
#        eval 'use Test::LongString';
#        skip('Test::LongString needed for these tests') if $@;
#        my $expected_abh = slurp('t/sample_data_abh.tab');
#        my $result_abh   = slurp('t/sample_data.tab_abh.tab');
#        is_string( $result_abh, $expected_abh, "ABH's corrected assigned" );
#    }
#}
