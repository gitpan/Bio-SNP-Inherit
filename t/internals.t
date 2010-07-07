use strict;
use warnings;
use Test::More tests=> 17;

BEGIN {
    use_ok('Bio::SNP::Inherit');
}

my @test_chars = qw{ T C };
my $expected_combined = 'CT';
my $result_combined = Bio::SNP::Inherit::_sort_and_join(@test_chars);
is($result_combined, $expected_combined, '_sort_and_join worked');

my $starting_word = 'howdy there';
my $expected_word = 'deehhortwy';
my $result_word = Bio::SNP::Inherit::_sorted_characters($starting_word);
is($result_word, $expected_word, '_sorted_characters worked');

my @expected_summary_data = qw{ 1  1  0  0  0  1  0  2  0  1  1                  }; 
my @result_summary_data = Bio::SNP::Inherit::_summarize_genotypes( qw{ AA AC GG -- GG TT CG } );
is_deeply(\@result_summary_data, \@expected_summary_data, 'summary line');


#This will not show up in the Devel::Cover report
eval{Bio::SNP::Inherit::_nonredundant_chars();};
like($@, qr/arguments required/, 'lack of arguments throws exception');

test_nonredundant_function( ['AAA'], ['A' ]);

test_nonredundant_function( ['AG'],  ['A', 'G'] );

test_nonredundant_function( [ 'AG', 'AG' ], ['A', 'G'] );
test_nonredundant_function( [ 'AA', 'AG', 'AA' ], ['A', 'G'] );

test_nonredundant_function( [ 'AGA', 'AA' ],      ['A', 'G'] );
test_nonredundant_function( [ '-A',  '-A' ],      ['A' ]);
test_nonredundant_function( [ 'A-',  'A-' ],      ['A' ]);
test_nonredundant_function( [ '--',  '--' ],      [q{-} ]);

test_nonredundant_function( [ 'CC', 'CC', 'CC' ], ['C' ]);
test_nonredundant_function( [ 'A',  'AA', 'AA' ], ['A' ]);

eval{Bio::SNP::Inherit::_nonredundant_chars('YZ','YZ');};
like($@,qr/Unknown character/, 'unknown character throws exception');

my @expected = qw{ A B C D E };
my @result = Bio::SNP::Inherit::_chars_from('ABCDE');

is_deeply(\@result, \@expected, '_chars_from');

sub test_nonredundant_function {
    my $aref     = shift;
    my $expected = shift;
    my @results = Bio::SNP::Inherit::_nonredundant_chars( @{$aref} );
    my $message = _quote_and_tab($aref) . " results in " . _quote_and_tab($expected);
    is_deeply(\@results, $expected, $message); 
}

sub _quote_and_tab {
    my $aref      = shift;
    my @elements  = @{$aref};
    my $QUOTE     = q{'};
    my $TAB_QUOTE = $QUOTE . "\t" . $QUOTE;
    return $QUOTE . join( $TAB_QUOTE, @elements ) . $QUOTE;
}
