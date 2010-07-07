
use strict;
use warnings;

use Module::Build 0.3601;


my %module_build_args = (
  'build_requires' => {
    'Module::Build' => '0.3601'
  },
  'configure_requires' => {
    'ExtUtils::MakeMaker' => '6.31',
    'Module::Build' => '0.3601'
  },
  'dist_abstract' => 'Module for determining the parental origin of specific SNPs based on genotype data.',
  'dist_author' => [
    'Christopher Bottoms <molecules at cpan.org>'
  ],
  'dist_name' => 'Bio-SNP-Inherit',
  'dist_version' => '0.00100001',
  'license' => 'perl',
  'module_name' => 'Bio::SNP::Inherit',
  'recommends' => {},
  'recursive_test_files' => 1,
  'requires' => {
    'Moose' => '0.90',
    'Test::More' => '0',
    'namespace::autoclean' => '0'
  },
  'script_files' => []
);


my $build = Module::Build->new(%module_build_args);

$build->create_build_script;
