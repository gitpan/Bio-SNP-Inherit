#!/usr/bin/perl
use strict;
use warnings;
use Test::More tests=> 9;

BEGIN{ use_ok('Bio::SNP::Inherit'); }
my $F1_gt = 'AA';
my $A_gt = ['AA', 'AA'];

my $B_gt = ['AA', 'AA'];

test_can_be_F1_of( 'AA', ['AA', 'AA'], ['AA', 'AA'], 1);
test_can_be_F1_of( 'AG', ['AG', 'AG'], ['AA', 'AA'], 1);
test_can_be_F1_of( 'AA', ['GG', 'GG'], ['AA', 'AA'], 0);
test_can_be_F1_of( 'AA', ['GG', 'GG'], ['CC', 'CC'], 0);
test_can_be_F1_of( 'AA', ['AG', 'AG'], ['GG', 'GG'],0);
test_can_be_F1_of( '--', ['AA', 'AA'], ['AA', 'AA'], 1);
test_can_be_F1_of( '--', ['TT', '--', 'TT'], ['--', '--', '--'], 1);
test_can_be_F1_of( 'TT', ['AA', 'AA', 'AA'], ['--', '--', '--'], 0);


sub test_can_be_F1_of {
    my $F1     = shift;
    my $A_aref = shift;
    my $B_aref = shift;
    my $expected = shift;

    my $operative_word = 'can';

    if($expected eq 0){
        $operative_word .= 'not';
    }

    my $result = Bio::SNP::Inherit::_can_be_F1_of($F1, $A_aref, $B_aref );
    my $message = "'$F1' $operative_word result from " 
                    . _space_and_tab($A_aref)
                    . ' and ' . _space_and_tab($B_aref);

    is($result, $expected, $message); 
}

sub _space_and_tab {
    my $aref      = shift;
    my @elements  = @{$aref};
    my $QUOTE     = q{'};
    my $TAB_QUOTE = $QUOTE . ' ' . $QUOTE;
    return $QUOTE . join( $TAB_QUOTE, @elements ) . $QUOTE;
}
