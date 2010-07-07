#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests=> 15;

BEGIN{ use_ok('Bio::SNP::Inherit'); }

can_ok('Bio::SNP::Inherit', '_calculate_abh_from_F1');

my %explanation = (
    'A' => 'parentA genotype',
    'B' => 'parentB genotype',
    'H' => 'heterozygous genotype',
    '-' => 'missing data',
    '--' => 'missing parental data',
    '!' => 'nonparental',
    '!!' => 'nonparental F1',
    '~' => 'nonpolymorphic',
    '#' => 'nonpolymorphic conflict',
    '%' => 'polymorphic parent',
);

test_remove_matching_char('C', ['A', 'C', 'G'], ['A', 'G']); 
test_remove_matching_char('A', ['A', 'C', 'G'], ['C', 'G']); 

sub test_remove_matching_char{
    my $char = shift;
    my $starting_aref = shift;
    my $expected_aref = shift;

    my @starting_array = @{ $starting_aref };
    my @result_array = Bio::SNP::Inherit::_remove_matching_char($char, @starting_array);
is_deeply(\@result_array, $expected_aref, '_remove_matching_char');
    
}



test_abh_from_F1('CG',['CC','CC'], ['CG'],'H');
test_abh_from_F1('CC',['CC','CC'], ['CG'],'A');
test_abh_from_F1('GG',['CC','CC'], ['CG'],'B');
test_abh_from_F1('GG',['GG','GG'], ['GG'],'~');
test_abh_from_F1('CC',['GG','GG'], ['GG'],'#');
test_abh_from_F1('CC',['GC','GG'], ['GG'],'%');
test_abh_from_F1('CC',['GG','GG'], ['GG', 'CC', 'AA'],'!!');
test_abh_from_F1('CC',['GG','GG'], ['GG', 'CC'],'!!');
test_abh_from_F1('CC',['GG','GG'], ['--', '--'],'--');
test_abh_from_F1('CC',['GC','GG'], ['GG'],'%');
test_abh_from_F1('AA',['GC','GG'], ['GG'],'!');

sub test_abh_from_F1 {
    my $sample_gt    = shift;
    my $parentA_aref = shift;
    my $F1_aref      = shift;
    my $expected     = shift;

    my $result = Bio::SNP::Inherit::_calculate_abh_from_F1( $sample_gt, $parentA_aref, $F1_aref );

    my $parentA_values = "'" . join( "'\t'", @{$parentA_aref} ) . "'";
    my $F1_values      = "'" . join( "'\t'", @{$F1_aref} ) . "'";
    my $why            = $explanation{$expected};

    is( $result, $expected,
              "\nSample genotype: '$sample_gt'     \n"
            . "ParentA values:   $parentA_values \n"
            . "ParentF1 values:  $F1_values \n"
            . "Expected result:  $expected ($why)\n"
            . "Actual result:    $result         \n" );

}
