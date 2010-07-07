use strict;
use warnings;
use Test::More tests => 6;

#Test that Bio::SNP::Inherit is loadable
BEGIN {
    use_ok( 'Bio::SNP::Inherit' );
}


#Test that Bio::SNP::Inherit has a 'manifest_filename' function/method
can_ok('Bio::SNP::Inherit','manifest_filename');

#Test that manifest_filename is a required attribute
eval { my $snp_obj = Bio::SNP::Inherit->new();};
like($@, qr/required/xms, "'manifest_filename' is required in the constructor");

my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/sample_manifest.tab',
data_filename => 't/sample_data.tab',);

isa_ok( $snp_obj, 'Bio::SNP::Inherit' );

_remove_temp_file('t/sample_data.tab_summary.tab');
_remove_temp_file('t/sample_data.tab_abh.tab');

sub _remove_temp_file{
    my $filename = shift;
    my $removed = unlink $filename;
    ok($removed, "removed tempoary file '$filename'");
}
