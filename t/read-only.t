use strict;
use warnings;
use Test::More tests => 6;

BEGIN {
    use_ok('Bio::SNP::Inherit');
}
my $snp_obj = Bio::SNP::Inherit->new(
    manifest_filename => 't/sample_manifest.tab',
    data_filename     => 't/sample_data.tab'
);
my @attribute_names = qw{
    _data_fh
    data_filename
    manifest_filename
    _sample_for
    sample_ids
};

for my $name (@attribute_names) {
    my $attribute = $snp_obj->meta->get_attribute($name);
    my $is_readable = ! ! $attribute->get_read_method();
    my $is_not_writeable = ! $attribute->get_write_method();
    my $is_read_only = $is_readable && $is_not_writeable;
    ok(  $is_read_only, "'$name' is read-only" );
}

