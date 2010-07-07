
use strict;
use warnings;



use ExtUtils::MakeMaker 6.31;



my %WriteMakefileArgs = (
  'ABSTRACT' => 'Module for determining the parental origin of specific SNPs based on genotype data.',
  'AUTHOR' => 'Christopher Bottoms <molecules at cpan.org>',
  'BUILD_REQUIRES' => {
    'Module::Build' => '0.3601'
  },
  'CONFIGURE_REQUIRES' => {
    'ExtUtils::MakeMaker' => '6.31',
    'Module::Build' => '0.3601'
  },
  'DISTNAME' => 'Bio-SNP-Inherit',
  'EXE_FILES' => [],
  'LICENSE' => 'perl',
  'NAME' => 'Bio::SNP::Inherit',
  'PREREQ_PM' => {
    'Moose' => '0.90',
    'Test::More' => '0',
    'namespace::autoclean' => '0'
  },
  'VERSION' => '0.00100001',
  'test' => {
    'TESTS' => 't/*.t'
  }
);


unless ( eval { ExtUtils::MakeMaker->VERSION(6.56) } ) {
  my $br = delete $WriteMakefileArgs{BUILD_REQUIRES};
  my $pp = $WriteMakefileArgs{PREREQ_PM};
  for my $mod ( keys %$br ) {
    if ( exists $pp->{$mod} ) {
      $pp->{$mod} = $br->{$mod} if $br->{$mod} > $pp->{$mod};
    }
    else {
      $pp->{$mod} = $br->{$mod};
    }
  }
}

delete $WriteMakefileArgs{CONFIGURE_REQUIRES}
  unless eval { ExtUtils::MakeMaker->VERSION(6.52) };

WriteMakefile(%WriteMakefileArgs);



