use strict;
use warnings;
use Test::More tests=> 27;

BEGIN{
use_ok('Bio::SNP::Inherit');
}


my %explanation = (
    'A' => 'parentA genotype',
    'B' => 'parentB genotype',
    'H' => 'heterozygous genotype',
    '-' => 'missing data',
    '--' => 'missing parental data',
    '!' => 'nonparental',
    '~' => 'nonpolymorphic',
    '#' => 'nonpolymorphic conflict',
    '%' => 'polymorphic parent',
);


#        sample     parentA           parentB               expected
test_abh( 'GG',     [ 'AA', 'AA', 'AA' ],     [ '--', '--', '--' ],     '--'     );
test_abh( 'AA',     [ 'AA', 'AA' ],     [ 'AA', 'AA' ],     '~'     );
test_abh( 'AA',     [ 'AA', 'AA' ],     [ '--', '--' ],     '--'     );
test_abh( 'AA',     [ '--', '--' ],     [ 'AA', 'AA' ],     '--'     );
test_abh( 'AA',     [ '--', '--' ],     [ '--', '--' ],     '--'     );
test_abh( 'CT',     [ 'CC', 'CC' ],     [ 'TT'       ],     'H'     );
test_abh( 'CC',     [ 'CT', 'CC' ],     [ 'TT'       ],     '%'     );
test_abh( 'CC',     [ 'TT', 'CC' ],     [ 'TT'       ],     '%'     );
test_abh( 'CC',     [ 'CT', 'CC', 'CC' ],     [ 'TT' ],     '%'     );
test_abh( 'CC',     [ 'CT', 'CC', 'CC' ],     [ 'CC' ],     '%'     );
test_abh( 'CC',     [ 'AG', 'AG', 'AG' ],     [ 'GG' ],     '!'     );
test_abh( 'AG',     [ 'AG', 'AG', 'AG' ],     [ 'AG' ],     '%'     );
test_abh( 'CC',     [ 'AG', 'AG', 'AG' ],     [ 'GG', 'GC' ],     '%'     );

#        sample     parentA   parentB       expected
test_abh( '--',     [ '--' ],   [ '--' ],   '-'     );
test_abh( 'AA',     [ 'A-' ],   [ 'AA' ],   '~'     );
test_abh( 'GC',     [ 'GG' ],   [ 'CC' ],   'H'     );
test_abh( 'GG',     [ 'GG' ],   [ 'CC' ],   'A'     );
test_abh( 'CC',     [ 'GG' ],   [ 'CC' ],   'B'     );
test_abh( 'GG',     [ 'GG' ],   [ 'TT' ],   'A'     );
test_abh( 'TG',     [ 'GG' ],   [ 'TT' ],   'H'     );
test_abh( 'TT',     [ 'GG' ],   [ 'TT' ],   'B'     );
test_abh( 'GG',     [ 'TT' ],   [ 'TT' ],   '#'     );
test_abh( 'GG',     [ 'GG' ],   [ 'GG' ],   '~'     );
test_abh( 'CC',     [ 'GG' ],   [ 'GG' ],   '#'     );
test_abh( 'AG',     [ 'GG' ],   [ 'TT' ],   '!'     );
test_abh( 'TA',     [ 'GG' ],   [ 'TT' ],   '!'     );

#AA AA AA (GG) -- -- --
sub test_abh{
    my $genotype = shift;
    my $parentA_aref = shift;
    my $parentB_aref = shift;
    my $expected = shift;

    my $result = Bio::SNP::Inherit::_calculate_abh(
                                                $genotype, 
                                                $parentA_aref,
                                                $parentB_aref
                                                );
    my $parentA_values = "'" . join("'\t'", @{ $parentA_aref }) . "'";
    my $parentB_values = "'" . join("'\t'", @{ $parentB_aref }) . "'";
    my $why = $explanation{$expected};
    is( $result,
        $expected, 
        "\nSample genotype: '$genotype'     \n"
        . "ParentA values:  $parentA_values \n"
        . "ParentB values:  $parentB_values \n"
        . "Expected result: $expected ($why)\n"
        . "Actual result:   $result         \n"
      );
}
