package Bio::SNP::Inherit;
BEGIN {
  $Bio::SNP::Inherit::VERSION = '0.00100001';
}
use strict;
#ABSTRACT: Module for determining the parental origin of specific SNPs based on genotype data.
use Moose;                  #Obect system. Also turns on strict and warnings.
use namespace::autoclean;   #"Keep imports out of your namespace"
use IO::File;               #Provides file methods (in lieu of <> )
use Carp qw{ confess };     #Die with stack trace
use List::MoreUtils qw{ all any none uniq zip }; #list-related functions
#use Smart::Comments '###'; #Used for debugging

#****************************************************************************#
# BUILD                                                               (START)#
# (This allows us to perform some tasks just after the object is created     #
#  but before it is returned to the user.)                                   #
#----------------------------------------------------------------------------#
sub BUILD{
    my $self = shift;

    #Get manifest data. Should die if there is an error.
    $self->_process_manifest_file();

    #Process data file and output summary and abh results
    $self->_output_summary_and_abh_files();

    return;
}
#----------------------------------------------------------------------------#
# BUILD                                                                 (END)#
#****************************************************************************#


#****************************************************************************#
# CONSTANTS                                                           (START)#
#----------------------------------------------------------------------------#
my $EMPTY_STRING = q{};

#Inheritance constants
my $MISSING           = q{-};
my $MISSING_PARENT    = $MISSING x 2;
my $NONPARENTAL       = q{!};
my $NONPARENTAL_F1    = $NONPARENTAL x 2;
my $POLYM             = q{%};
my $ALLELE_A          = q{A};
my $ALLELE_B          = q{B};
my $HETEROZYGOUS      = q{H};
my $NONPOLYM          = q{~};
my $NONPOLYM_CONFLICT = q{#};
#----------------------------------------------------------------------------#
# CONSTANTS                                                             (END)#
#****************************************************************************#


#****************************************************************************#
# MANIFEST FILE                                                       (START)#
#----------------------------------------------------------------------------#
has 'manifest_filename' => (
    is=>'ro',
    isa=>'Str',
    required => 1,
);

has '_sample_for' => (
    traits => ['Hash'],
    is      => 'ro',
    isa     => 'HashRef',
    default => sub{ {} },
    handles =>{
        _get__sample_for => 'get',
        _sample_ids_from_manifest => 'keys',
    },
);

sub _process_manifest_file {
        my $self = shift;
        my %sample_for;

        #Open filehandle from manifest filename
        my $fh = IO::File->new($self->manifest_filename, '<') or die $@;

        #discard first line
        my $first_line = <$fh>;

        #read in each id with its associated data
        while ( my $line = <$fh> ) {

            #remove newline character
            $line = _sans_newlines($line);

            #load each value from the line into appropriate variables
            my ( $id, $name, $group, $parentA, $parentB, $replicate_of,
                $F1_ancestor )
                = split /\t/, $line;

            $id = _trim($id);

            #confess or next . . .
            #If $id is an empty string, throw an error if anything but
            #   whitespace exists on the line
            if ( $id eq $EMPTY_STRING ) {
                if ( $line =~ m{ \A \s+ \z}xms ) {
                    #Skip this empty line
                    next;
                }
                else {
                    confess "Empty id in manifest file: '$line'";
                }
            }

            #confess if hash for this given id already exists
            confess 'duplicate sample id in manifest file' . "'$id'"
                if exists $sample_for{$id};

            #remove leading/trailing whitespace
            foreach ( $id, $name, $group, $parentA, $parentB, $replicate_of, $F1_ancestor ) {
                $_ = _trim($_);
            }

            $sample_for{$id} = {};
            my $sample = $sample_for{$id};

            #only store values that aren't empty
            _add_if_not_empty( $sample, name         => $name );
            _add_if_not_empty( $sample, group        => $group );
            _add_if_not_empty( $sample, parentA       => $parentA );
            _add_if_not_empty( $sample, parentB       => $parentB );
            _add_if_not_empty( $sample, replicate_of => $replicate_of );
            _add_if_not_empty( $sample, F1_ancestor => $F1_ancestor );

            #If this is a replicate of another sample, then store id in
            # replicates array of the proper sample
            if($replicate_of){
                push @{ $sample_for{$replicate_of}{replicates} }, $id;
            }
        }

        #Finished reading data, so close the filehandle
        $fh->close();

        #Store reference to the hash containing all of the data
        $self->{_sample_for} = \%sample_for;

        return;
}

#Trim away leading and trailing whitespace
sub _trim{
    my $string = shift;

    #Don't try substitutions on undef
    return if ! defined( $string);

    #Trim all leading whitespace
    $string =~ s/\A \s+//xms;

    #Trim all trailing whitespace
    $string =~ s/\s+ \z//xms;

    #Return trimmed string
    return $string;
}

#Adds key, value pairs to a hash if the value is not empty.
sub _add_if_not_empty{
    my $hashref = shift;
    my $key = shift;
    my $value = shift;

    #Only assign to hash for keys that are defined and not empty
    if(defined($value) && $value ne $EMPTY_STRING){
        ${$hashref}{$key} = $value;
    }
    return $hashref;
}
#----------------------------------------------------------------------------#
# MANIFEST FILE                                                         (END)#
#****************************************************************************#


#****************************************************************************#
# SUMMARIZE DATA                                                      (START)#
#----------------------------------------------------------------------------#

#Must be specified in the constructor. This is indicated here by (1) making it
# required and by (2) not giving it a default nor a builder
has 'data_filename' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
);

# leading underscore indicates a private attribute by convention
has '_data_fh' => (
    is       => 'ro',
    isa      => 'IO::File',
    lazy     => 1,
    required => 1,
    builder  => '_build__data_fh',
);

sub _get_data_line{
    my $self = shift;
    return $self->_data_fh->getline();
}

#sub name = '_build_' . '_data_fh'
sub _build__data_fh {
    my $self = shift;
    my $fh = IO::File->new( $self->data_filename(), '<' ) or die $@;
    return $fh;
}

has 'sample_ids' => (
    traits   => ['Array'],
    is       => 'ro',
    isa      => 'ArrayRef',
    lazy     => 1,
    required => 1,
    builder  => '_build_sample_ids',
    handles  => {
        _ids_from_data_file => 'elements',
    },
);

sub _build_sample_ids {
    my $self = shift;

    #skip all lines until finding the [Data] label.
    while ( my $line = $self->_data_fh->getline() ) {
        last if $line =~ /\[ data \]/xmsi;    #i makes it case insensitive
    }

    #read in line containing ids
    my $line_of_ids = $self->_data_fh->getline();
    if ( !defined $line_of_ids ) {
        my $dying_message =
            'either no sample ids or no "[Data]" tag found in file: '
            . $self->data_filename;
        confess $dying_message;
    }

    #remove newline(s)
    $line_of_ids = _sans_newlines($line_of_ids);

    #extract ids from line, ignoring blank field at the beginning
    my ( undef, @ids ) = split /\t/, $line_of_ids;
    return \@ids;
}

has '_columns_for' => (
    traits  => ['Hash'],
    is      => 'ro',
    isa     => 'HashRef',
    lazy    => 1,
    builder => '_build__columns_for',
    handles => { _columns_for_id_aref => 'get', },
);

#getting columns corresponding to given id.
sub _get_columns_for{
    my $self = shift;
    my $id = shift;
    confess if ! defined($id);
    return @{ $self->_columns_for_id_aref($id) };
}

sub _build__columns_for{
    my $self = shift;

    #get all of the sample ids, in their column order
    my @data_ids = $self->_ids_from_data_file();

    #get all of the column indices
    my @indices = (0 .. $#data_ids);

    #store column indicdes
    my @columns_aref;
    for my $index (0 .. $#indices){
        push @columns_aref, [ $index ];
    }

    #Create hash to hold column info
    my %columns_for = zip(@data_ids, @columns_aref);

    #-------------------------------------------------------------------------
    #For each sample that has replicates, add column indices of the replicates
    #  This just checks ids from the manifest file, since they are the only
    #  ones that have replicate information.
    for my $manifest_id ( $self->_sample_ids_from_manifest() ) {

        my $sample = $self->{_sample_for}->{$manifest_id};

        #if sample has replicates, add column indices of replicates
        if ( defined $sample->{replicates} ) {
            for my $replicate ( @{ $sample->{replicates} } ) {
                if( defined $columns_for{$replicate}){
                push @{ $columns_for{$manifest_id} }, @{ $columns_for{$replicate} };
                }
            }
        }
    }
    #-------------------------------------------------------------------------

    return \%columns_for;
}

sub _in_group{
    my $self = shift;
    my $group_regex = shift;
    my $id = shift;
    my $sample = $self->{_sample_for}->{$id};
    return 0 if(!defined $sample);
    if($sample->{group} =~ / $group_regex /xms){
        return 1;
    }else{
        return 0;
    }
}

sub _parented_href{
    my $self = shift;
    my %parented = (
        names => [],
        ids => [],
        groups => [],
        parentA_names => [],
        parentB_names => [],
        F1_ancestors  => [],
    );


    GROUPS_OF_CODES:
    for my $id ($self->_ids_from_data_file()) {
        my $sample = $self->{_sample_for}->{$id};

        #Skip this id, if no parentA
        next if ! defined( $sample->{parentA});

        if ( defined( $sample->{parentB} ) || defined( $sample->{F1_ancestor} )) {

            my $parentA_id = $sample->{parentA};
            my $parentA_name = $self->{_sample_for}->{$parentA_id}->{name};

            #Get $parentB_id.
            #Change it to undef if it is the empty string
            my $parentB_id = $sample->{parentB} || undef;

            my $F1_ancestor_id = $sample->{F1_ancestor} || undef;

            my $parentB_name;
            if ($parentB_id) {
                $parentB_name = $self->{_sample_for}->{$parentB_id}->{name};
            }
            elsif ($F1_ancestor_id) {
                $parentB_name =
                    '(F1)' . $self->{_sample_for}->{$F1_ancestor_id}->{name};
            }

            push @{ $parented{ids} },           $id;
            push @{ $parented{parentA_names} }, $parentA_name;
            push @{ $parented{parentB_names} }, $parentB_name;
            push @{ $parented{groups} },        $sample->{group};
            push @{ $parented{names} },         $sample->{name};
        }
        else{
            #Skip this one since it doesn't have a parentB or an F1 ancestor
            next;
        }
    }
    return \%parented;
}

sub _output_summary_and_abh_files {
    my $self = shift;

    #Column headers for Summary output file
    $self->_summary_fh->print($self->_summary_header());

    my @all_sample_ids = $self->_ids_from_data_file;

    my $nam_regex = qr{ NAM [ ]* F1}xms;
    my @nam_f1_ids =
        grep { $self->_in_group( $nam_regex, $_ ) } @all_sample_ids;

    my %parented = %{ $self->_parented_href() };

    #Columns headers for ABH output file
    for my $item (qw{ ids parentA_names parentB_names groups names}) {
        my $abh_header_row = "\t" . join "\t", @{ $parented{$item} };
        $self->_abh_fh->print($abh_header_row . "\n");
    }

    #Read data from file and analyze it one line at a time
    # (if already done, this will be skipped due to filehandle
    # already being at the end of the file)
    while ( my $line = $self->_get_data_line() ) {

        #Remove newline character
        $line = _sans_newlines($line);

        #Extract SNP name and the data from this line
        my ( $snp, @genotypes ) = split /\t/, $line;

        #Remove whitespace at beginning and end of SNP name
        $snp = _trim($snp);

        #Skip blank lines
        next if $snp eq $EMPTY_STRING;

        #Print SNP name in first column of current line of ABH output file
        $self->_abh_fh->print($snp);

        #Tally occurrences of eacy genotype
        my @summarized_genotypes = _summarize_genotypes(@genotypes);

        #Associated data ids with genotypes
        my %genotype_for = zip( @all_sample_ids, @genotypes );

        #Initialize counters with zeros. Thus zeros will show up in reports.
        my $inconsistent_count = 0;
        my %num_of = (
            parentA => {
                $POLYM   => 0,
                $MISSING => 0,
            },
            parentB => {
                $POLYM   => 0,
                $MISSING => 0,
            },
            F1 => { $MISSING => 0, },
        );

        #----------------------------------------------------------------START
        #Summarize this one line of data for all of the NAM F1's
        #---------------------------------------------------------------------
    NAM_F1_LOOP:
        for my $sample_id (@nam_f1_ids) {

            my $sample    = $self->_get__sample_for($sample_id);
            my $sample_gt = $genotype_for{$sample_id};

           #For our current purposes, we will not evaluate an F1 unless there is data for both parents.
            my $both_parents_present =
                defined $sample->{parentA} && defined $sample->{parentB};

            next NAM_F1_LOOP if (! $both_parents_present);

            my @genotype_chars;
            my $parents_missing_data_for_this_sample_count;
            my %gt_chars_of;
            for my $parent (qw{ parentA parentB}) {
                my @parent_columns =
                    $self->_get_columns_for( $sample->{$parent} );
                my @parent_genotypes = @genotypes[@parent_columns];

                #_nonredundant_chars
                my @parent_gt_chars = @{ $gt_chars_of{$parent} } =
                    _nonredundant_chars(@parent_genotypes);
                push @genotype_chars, @parent_gt_chars;

                my $num_parent_char = @parent_gt_chars;
                if ( $num_parent_char > 1 ) { $num_of{$parent}{$POLYM}++; }
                elsif ( $num_parent_char == 1 ) {
                    if ( $parent_gt_chars[0] eq $MISSING ) {
                        $num_of{$parent}{$MISSING}++;
                        $parents_missing_data_for_this_sample_count++;
                    }
                }
                elsif ( $_ < 1 ) {
                    confess
                        '_nonredundant_chars should return at least one char';
                }
            }

            #Keep track of how many F1s are missing data
            if ( $sample_gt eq $MISSING x 2 ) {
                $num_of{F1}{$MISSING}++;
            }
            elsif (
                !_can_be_F1_of(
                    $sample_gt, $gt_chars_of{parentA},
                    $gt_chars_of{parentB}
                )
                )
            {
                $inconsistent_count++;
            }
        }
        my $sum_of_odds_and_ends = $inconsistent_count
                                    + $num_of{parentA}{$POLYM}
                                    + $num_of{parentB}{$POLYM}
                                    + $num_of{parentA}{$MISSING}
                                    + $num_of{parentB}{$MISSING};

        my $summary_line = join "\t", $snp, @summarized_genotypes,
            $inconsistent_count, $num_of{parentA}{$POLYM},
            $num_of{parentB}{$POLYM},   $num_of{parentA}{$MISSING},
            $num_of{parentB}{$MISSING}, $sum_of_odds_and_ends,
            $num_of{F1}{$MISSING};

        $self->_summary_fh->print($summary_line . "\n");
        #---------------------------------------------------------------------
        #Summarize this one line of data for all of the NAM F1's
        #------------------------------------------------------------------END


        #Output inheritance codes
        #First step: calculate for those with two parents
        #Second step: calculate for those with parentA and one F1 ancestor
        my $abh_line = $EMPTY_STRING;

        #----------------------------------------------------------------START
        #Calculate inheritance codes for everything with two parents
        #---------------------------------------------------------------------
        for my $sample_id ( @{ $parented{ids} } ) {

            my $sample    = $self->_get__sample_for($sample_id);
            my $sample_gt = $genotype_for{$sample_id};
            my $parentA_id        = $sample->{parentA};
            my @parentA_columns   = $self->_get_columns_for($parentA_id);
            my @parentA_genotypes = @genotypes[@parentA_columns];

            my $sample_abh;
            if ( $sample->{parentB} ) {
                my $parentB_id        = $sample->{parentB};
                my @parentB_columns   = $self->_get_columns_for($parentB_id);
                my @parentB_genotypes = @genotypes[@parentB_columns];

                $sample_abh = _calculate_abh( $sample_gt, \@parentA_genotypes,
                    \@parentB_genotypes );
            }elsif($sample->{F1_ancestor}){
                my $F1_ancestor_id        = $sample->{F1_ancestor};
                my @F1_columns   = $self->_get_columns_for($F1_ancestor_id);
                my @F1_genotypes = @genotypes[@F1_columns];
                $sample_abh = _calculate_abh_from_F1( $sample_gt, \@parentA_genotypes,
                    \@F1_genotypes);
            }
            $abh_line .= "\t" . $sample_abh;
        }

        #print abh_line to file and add newline
        $self->_abh_fh->print($abh_line ."\n");
        #------------------------------------------------------------------END


    }
    $self->_summary_fh->flush();
    $self->_summary_fh->close();
    $self->_abh_fh->flush();
    $self->_abh_fh->close();

    return;
}


has '_abh_fh' => (
    is       => 'ro',
    isa      => 'IO::File',
    lazy     => 1,
    required => 1,
    builder  => '_build__abh_fh',
);

sub _build__abh_fh {
    my $self     = shift;
    my $filename = $self->data_filename() . '_abh.tab';
    my $fh       = IO::File->new( $filename, '>' ) or die $@;
    return $fh;
}

has '_summary_fh' => (
    is       => 'ro',
    isa      => 'IO::File',
    lazy     => 1,
    required => 1,
    builder  => '_build__summary_fh',
);

sub _build__summary_fh {
    my $self     = shift;
    my $filename = $self->data_filename() . '_summary.tab';
    my $fh       = IO::File->new( $filename, '>' ) or die $@;
    return $fh;
}

sub _sorted_characters{
    my $string = shift;
    $string =~ s/ \s //xms;
    #mea culpa
   return join $EMPTY_STRING, sort split $EMPTY_STRING, $string;

}

sub _sort_and_join{
    my @chars = @_;
    return join $EMPTY_STRING, sort @chars;
}

sub _can_be_F1_of{
    my $sample_gt = shift;
    my $parentA_aref = shift;    # 'aref' means array reference
    my $parentB_aref = shift;
    my @parentA_chars = _nonredundant_chars( @{$parentA_aref} );
    my @parentB_chars = _nonredundant_chars( @{$parentB_aref} );

    #If parentA or parentB only have missing data, then they could be
    #   anything
    if (all { $_ eq $MISSING } @parentA_chars){
        @parentA_chars = qw{ A C G T - };
    }
    if (all { $_ eq $MISSING } @parentB_chars){
        @parentB_chars = qw{ A C G T - };
    }

    my $sorted_sample_gt = _sorted_characters($sample_gt);
    return 1 if $sorted_sample_gt eq $MISSING x 2;

    for my $parentA_char (@parentA_chars){
       for my $parentB_char (@parentB_chars){
            my $expected_gt = _sort_and_join($parentA_char, $parentB_char);
            return 1 if($expected_gt eq $sorted_sample_gt);
       }
    }

    #If nothing matched, return false
    return 0;
}


sub _calculate_abh_from_F1 {
    my $sample_gt    = shift;
    my $parentA_aref = shift;
    my $F1_aref      = shift;

    _calculate_abh($sample_gt, $parentA_aref, $F1_aref, 'F1');
}

sub _remove_matching_char{
    my $char = shift;
    my @array = @_;
    my @final_array;
    for(@array){
        next if $_ eq $char;
        push @final_array, $_;
    }
    return @final_array;
}

sub _equivalent_arrays{
    my $first_aref = shift;
    my $second_aref = shift;
    my @first = @{ $first_aref };
    my @second = @{ $second_aref };

    my $same_size = scalar @first == scalar @second;
    return if ! $same_size;

    for (0 .. $#first ){
        return if $first[$_] ne $second[$_] ;
    }
    return 1;
}

sub _calculate_abh {
    my $sample_gt     = shift;    # 'gt' means genotype
    my $parentA_aref  = shift;    # 'aref' means array reference
    my $parentB_aref  = shift;
    my $parentB_is_F1 = shift;
    my @parentA_chars = _nonredundant_chars( @{$parentA_aref} );
    my @parentB_chars = _nonredundant_chars( @{$parentB_aref} );

    #If genotype has one or more $MISSING alleles, then return $MISSING
    return $MISSING if ( $sample_gt =~ m{ $MISSING }xms );

    #Find out if the parents have exactly the same genotype
    my $parents_same;
    if ( _equivalent_arrays(\@parentA_chars, \@parentB_chars )) {
        $parents_same = 1;
    }

    #Return 'missing' if one of the parents have only missing data
    #Return error code if sample has alleles not found in parents
    #Must check in this order, otherwise having no parental info will be
    #   misdiagnosed as having a nonparental allele
    if ( @parentA_chars == 1 && $parentA_chars[0] eq $MISSING ) {
        return $MISSING_PARENT;
    }
    elsif ( @parentB_chars == 1 && $parentB_chars[0] eq $MISSING ) {
        return $MISSING_PARENT;
    }
    else {

        #combine parental character sets
        my @both_parent_chars =
            _nonredundant_chars( @parentA_chars, @parentB_chars );

        #Check if the sample genotype can have descended from parents
        if ( !_is_comprised_from( $sample_gt, @both_parent_chars ) ) {
            if ($parents_same) {
                return $NONPOLYM_CONFLICT;
            }
            else {
                return $NONPARENTAL;
            }
        }
        else {

            #Good! No errors and no missing data.
        }
    }

    if ( !$parentB_is_F1 ) {

       #If genotype only contains alleles from parents, then assign A, B, or H
        my $allele_A_count = 0;
        my $allele_B_count = 0;

        #Check if both parents are nonpolymorphic
        if ( @parentA_chars == 1 && @parentB_chars == 1 ) {
            my $parentA_char = $parentA_chars[0];
            my $parentB_char = $parentB_chars[0];

            #Indicate noninformative if parents are the same
            return $NONPOLYM if ($parents_same);

            #Extract characters from sample genotype
            my @sample_chars = _chars_from($sample_gt);

            #Tally how many are from which parent
            for my $sample_char (@sample_chars) {
                if($sample_char eq $parentA_char){
                    $allele_A_count++;
                    }
                elsif ($sample_char eq $parentB_char) { $allele_B_count++;}
            }

            if ( $allele_A_count == 2 ) {
                return $ALLELE_A;
            }
            elsif ( $allele_B_count == 2 ) {
                return $ALLELE_B;
            }
            elsif ( $allele_A_count == 1 && $allele_B_count == 1 ) {
                return $HETEROZYGOUS;
            }
            else {

                #Since we already checked for values that were missing or
                #   nonparental, this should be unreachable
                my $dying_message = "error calculating 'ABH' from:\n"
                    . join "'\n",
                    q{parentA => '} . join( "' '", @{$parentA_aref} ),
                    q{parentB => '} . join( "' '", @{$parentB_aref} ),
                    q{sample genotype =>'} . $sample_gt;

                confess $dying_message;
            }
        }
        else {

           #indicate that it is noninformative because of a polymorphic parent
            return $POLYM;
        }
    }
    else {
        for my $F1_gt ( @{$parentB_aref} ) {
            return $NONPARENTAL_F1
                if !_can_be_F1_of( $F1_gt, $parentA_aref, [$MISSING] );
        }

        return $POLYM if @parentA_chars > 1;

            if (@parentB_chars == 2) {
                my $parentA_char = $parentA_chars[0];
                my @revised_parentB_chars =
                    _remove_matching_char( $parentA_char, @parentB_chars );
                return _calculate_abh( $sample_gt, $parentA_aref,
                    \@revised_parentB_chars );
            }
            elsif (@parentB_chars == 1) {
                if ( _equivalent_arrays(\@parentB_chars, \@parentA_chars )) {
                    return _calculate_abh( $sample_gt, $parentA_aref,
                        $parentB_aref );
                }
                else {
                    ### $sample_gt
                    ### @parentA_chars
                    ### @parentB_chars
                    confess "I didn't think this was reachable";
                }
            }
    }
}

sub _summary_header{
    return "\t" . join("\t", qw{ AA AC AG AT CC CG CT GG GT TT --
                                inconsistent parentA_polymorphic
                                parentB_polymorphic parentA_unknown
                                parentB_unknown sum_of_odds_and_ends
                                F1_missing }) . "\n";
}

sub _summarize_genotypes {
    my @data      = @_;
    my @genotypes = qw{ AA AC AG AT CC CG CT GG GT TT -- };
    my %count;

    #Initialize hash %count, so that zeros show up in the summary
    for my $genotype (@genotypes) {
        $count{$genotype} = 0;
    }

    #Count occurrences of each genotype
    for my $datum (@data){

        #Remove any leading or trailing whitespace
        $datum = _trim($datum);

        #Get just the first two characters of each data string
        my $sorted_char_datum = _sorted_first_two_char($datum);

        #Add one to the count for this genotype
        $count{$sorted_char_datum}++;
    }

    #Return an array of just the specified genotypes
    #   (this is a 'hash slice')
    return @count{ @genotypes };
}
#----------------------------------------------------------------------------#
# SUMMARIZE DATA                                                        (END)#
#****************************************************************************#

#****************************************************************************#
# GENERAL UTILITIES                                                          #
#----------------------------------------------------------------------------#

sub _nonredundant_chars {
    my @strings            = @_;
    confess 'arguments required' if( @strings < 1);
    for( 0 .. $#strings){
        $strings[$_] = _sans_newlines($strings[$_]);
    }

    #Get just the unique characters out of the strings
    my $combined = join $EMPTY_STRING, @strings;
    my @single_chars = uniq _chars_from($combined);

    #Remove 'missing data' and empty characters
    my @filtered_chars =
        grep { $_ ne $MISSING && $_ ne $EMPTY_STRING } @single_chars;

    #Return $MISSING if there aren't any left
    return $MISSING if(@filtered_chars == 0);

    #Check that the characters are all DNA bases
    my $DNA_BASES = qr{ [ACGT] }xms;
    for my $filtered_char(@filtered_chars){

        #Throw exception if character is not a valid DNA base
        if($filtered_char !~ m{\A $DNA_BASES \z}xms){
            confess "Unknown character: '$filtered_char'";
        }
    }

    #Return nonredundant characters, now that they've all been checked
    return @filtered_chars;
}

#   Are all of the characters in a string a subset of the characters found in
#       this array of strings?
sub _is_comprised_from{
    #Get the string to be analyzed, followed by strings with acceptable chars
    my ($string, @acceptable_strings) = @_;

    #Put all of the characters from all of these strings into one array
    my @acceptable_chars;
    for my $acceptable_string( @acceptable_strings){
        push @acceptable_chars, uniq _chars_from($acceptable_string);
    }

    confess 'Empty string found'
         if _is_an_element_of($EMPTY_STRING,  @acceptable_chars);

    #Get all of the characters out of the string to be analyzed
    my @chars = uniq _chars_from($string);

    #If any character is not acceptable, return false
    for my $char (@chars){
       if(! _is_an_element_of($char, @acceptable_chars)){
            return 0;
       }
    }

    #Since no unacceptable character was found, return true
    return 1;
}

sub _is_an_element_of{

    my $query = shift;
    my @elements = @_;
    for my $element( @elements){
        return 1 if $query eq  $element;
    }

    return 0;
}

sub _chars_from{
    my @strings = shift;

    my @string_chars;
    for my $string (@strings){
        push @string_chars, split $EMPTY_STRING, $string;
    }

    #Return the individual characters from the string as an array
    return @string_chars;
}

sub _sorted_first_two_char {
    my $string = shift;

    #Get first two characters from the string
    my $first_two_char = substr $string, 0, 2;

    #Sort these two characters alphabetically
    my @sorted_char = sort { $a cmp $b } split $EMPTY_STRING, $first_two_char;

    #Concatenate and return the characters
    return join $EMPTY_STRING, @sorted_char;
}

sub _sans_newlines{
    my $string = shift;
    chomp $string;

    $string =~ s{ \r }{}xms;

    my $string_sans_newlines = $string;

    return $string_sans_newlines;
}


__PACKAGE__->meta->make_immutable;

1; # End of Bio::SNP::Inherit

__END__
=head1 NAME

=head1 SYNOPSIS

=head1 VERSION

version 0.00100001

    my $foo = Bio::SNP::Inherit->new(
        manifest_filename => 'manifest.tab',
        data_filename     => 'data.tab'
    );

    #Upon object construction, this outputs a summary file
    #   'data.tab_summary.tab' and a detailed file 'data.tab_abh.tab'
    #   containing parental allele designations for each sample that has
    #   parents defined for it in the manifest file

=head1 DESCRIPTION

This is a module for converting Single Nucleotide Polymorphism (SNP) genotype
data to parental allele designations. This helps with creating files suitable
for mapping, identifying and characterizing crossovers, and also helps with
quality control.

=head1 SUBROUTINES/METHODS

=head2 BUILD

    Since the integrity of the data in the manifest file is absolutely vital,
    building an object fails if there are duplicate sample ids in the
    manifest file.

=head1 ATTRIBUTES

=head2 manifest_filename

    Name of the file containing information for each sample id

    Required in the constructor

    The first line contains headers and the remaining lines contain
        tab-delimited fields in the following order:

        sample id     or "Institute Sample Label"    (e.g. "WG0096796-DNAA05" )
        sample name   or "Sample name"               (e.g. "B73xB97"          )
        group name    or "Group"                     (e.g. "NAM F1"           )
        parentA       or "Mother"                    (e.g. "WG0096795-DNAA01" )
        parentB       or "Father"                    (e.g. "WG0096796-DNAF01" )
        replicate of  or "Replicate(s)"    (id of sample that this replicates 
                                              e.g. "WG0096796-DNAA05"         ) 
        AxB F1        or "F1 of parentA and parentB" (e.g. "WG0096795-DNAA02" )

    The last four fields can be blank, if they are not applicable. However,
        being blank when they are applicable will result in failure of the
        program to analyze the data properly

=head2 data_filename

    Name of the tab-delimited file containing the data to be processed.

    Required in the constructor.

    The text '[Data]' in a line indicates that remaining lines are all data.
    The next line contains column headers, which are in fact the sample ids.
        Sample ids missing from the manifest file will not be processed.
    The next line contains the name of the SNP in the first field and data in
        the remaining fields.

    Data must be in the format of SNP_name{tab}AA{tab}GG{tab}.

=head1 OUTPUT FILES

    Upon object construction, two files are produced: one that summarizes the
    input and another that that describes the genotypes of samples in terms of
    their "parents". For example, a sample with a genotype of "CG" whose
    'parentA' has a genotype of "CC" and whose 'parentB' has a genotype of
    "GG" would have a heterozygous genotype, labeled as 'H'.

    Here are the possible allele designations that result:

        Allele designations for informative genotypes:
            A = parentA genotype
            B = parentB genotype
            H = heterozygous genotype

        Allele designations for noninformative genotypes:
            ~ = nonpolymorphic parents (i.e. both parents have same genotype)
            - = missing data
            -- = missing data for at least one parental
            % = polymorphic parent

        Error codes:
            # = conflict of nonpolymorphic expectation, meaning both parents
                    have the same genotype, but the sample has a different
                    genotype. For example, parentA and parentB both have the
                    genotype 'CC', but the sample has a genotype of 'TT'.

            ! = nonparental genotype, meaning each parent has a different
                    genotype, but the sample has at least one allele not seen
                    in either parent. For example, getting 'AG' for the
                    offspring when the parents have 'GG' and 'TT'.
                    (This should not even be seen when the data was obtained
                    from a biallelic assay.)

            !! = genotype of the F1 for parentA x parentB is incongruent with
                    the genotype for parentA

    See the bundled tests for examples.

=head1 TODO

    Output report detailing which samples have been processed and in what way.
    Also give descendents and ancestor relationships.

    Document ability to process files using F1 and parentA info (i.e. in the
    absence of parentB info).

    Add simple means of adding map info so that distances and chromosomes are
    output along with the marker names.

    Give crossover info?

    Give introgressions/regions attributable to specific ancestor(s).

    Use benchmarking to find out which (if any) to memoize:
    _nonredundant_chars
    _trim
    _is_comprised_from
    _sorted_characters
    _sort_and_join
    _chars_from
    _sorted_first_two_char

    Test bad file names

=head1 DIAGNOSTICS

    TODO

=head1 CONFIGURATION AND ENVIRONMENT

   TODO

=head1 DEPENDENCIES

   TODO

=head1 INCOMPATIBILITIES

   TODO

=head1 BUGS

Please report any you find. None have been reported as of the current release.

=head1 LIMITATIONS

This is ALPHA code. Use at your own risk. There are some major changes
that I want to do to it.

Be consciencious with the preparation of your input files (i.e. manifest file
and data file). Correct results depend on correct input files.

=head1 AUTHOR

Christopher Bottoms, C<< <molecules at cpan.org> >>

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::SNP::Inherit

=head1 ACKNOWLEDGEMENTS

=head1 LICENSE AND COPYRIGHT

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

Copyright 2010 Christopher Bottoms.

=cut