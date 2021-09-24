#!/usr/bin/env perl
package MyApp {
    use MooseX::App qw(Color);
    use Log::Any '$log';

    has 'log' => (
        is            => 'ro',
        isa           => 'Log::Any::Proxy',
        required      => 1,
        default       => sub { Log::Any->get_logger },
        documentation => 'Keep Log::Any::App object',
    );

    __PACKAGE__->meta->make_immutable;
}

package MyApp::GetProperPairsAndOverlappingReads {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use namespace::autoclean;
    use Data::Printer;
    use List::MoreUtils qw( zip );

    has_file 'bed6' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required => 1,
        must_exist => 1,
        documentation => q[bed6 of all aligned reads to the human genome ],
    );

    has_file 'sam_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required => 1,
        must_exist => 1,
        documentation => q[sam file without the header with all aligned reads],
    );


    has 'insertmax' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required => 1,
        is => 'ro',
        isa => 'Num',
        documentation => q[maximum distance between the parid end reads],
    );


    sub bed6_to_hash {

        my ( $self, $hash_ref ) = @_;
        my @header = qw/chr start end read_id score strand/;

        open( my $fh, '<', $self->bed6 ) or die "Could not open file '$self->bed6' $!";
        while ( my $row = <$fh> ) {
            chomp $row;
            my @split_row = split( "\t", $row );

            my $read_id;
            my $read_type;
            if ( $split_row[ 3 ] =~ m/(\S+)_(\S+)\/.*/ ) {
                $read_id   = $1;
                $read_type = $2;
            }

            my %bed6_row_hash = zip @header, @split_row;

            $hash_ref->{ $read_id }{ $read_type } = \%bed6_row_hash;

        }

    }


    sub sam_to_hash {

        my ( $self, $sam_file ) = @_;
        my %sam_hash;
        my @header = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL NUMOFDIF ALIGNSCORE/;

        open( my $fh, "<", $sam_file );
        while ( my $row = <$fh> ) {
            chomp $row;
            my @splitted_row = split( "\t", $row );

            my %row_hash;
            @row_hash{ @header } = @splitted_row;

            if ( $splitted_row[ 0 ] =~ m/^(\S+)_\S+/ ) {
                my $read_id = $1;
                $sam_hash{ $read_id }{ $splitted_row[ 0 ] } = \%row_hash;
                $sam_hash{ $read_id }{ selected } = "FALSE";
            }

        }

        return ( %sam_hash );

    }


    sub clean_single_read {
        my ( $self, $bed6_hash_ref ) = @_;

        foreach my $read_id ( sort { $a cmp $b } keys %{ $bed6_hash_ref } ) {

            my $length_pairs;
            $length_pairs = scalar keys %{ $bed6_hash_ref->{ $read_id } };

            if ( $length_pairs == 1 ) {
                delete $bed6_hash_ref->{ $read_id };
            }

        }

    }


    sub filter_reads {
        my ( $self, $sam_hash_ref, $bed6_hash ) = @_;

        foreach my $pair ( sort { $a cmp $b } keys %{ $bed6_hash } ) {

            if ( $bed6_hash->{ $pair }{ 'bait' }{ 'strand' } ne $bed6_hash->{ $pair }{ 'target' }{ 'strand' } ) {

                if ( $bed6_hash->{ $pair }{ 'bait' }{ 'strand' } eq "+" ) {

                    my $target_end = $bed6_hash->{ $pair }{ 'target' }{ 'end' };
                    my $bait_start = $bed6_hash->{ $pair }{ 'bait' }{ 'start' };

                    if ( $target_end  >= $bait_start and ( $target_end - $bait_start ) >= 15 and ( $target_end - $bait_start ) <= $self->insertmax  ) {
                        $sam_hash_ref->{ $pair }{ selected } = "TRUE";
                    }

                }
                elsif ( $bed6_hash->{ $pair }{ 'bait' }{ 'strand' } eq "-" ) {

                    my $bait_end = $bed6_hash->{ $pair }{ 'bait' }{ 'end' };
                    my $target_start = $bed6_hash->{ $pair }{ 'target' }{ 'start' };

                    if ( $bait_end >= $target_start and ( $bait_end - $target_start  ) >= 15 and ( $bait_end - $target_start ) <= $self->insertmax) {
                        $sam_hash_ref->{ $pair }{ selected } = "TRUE";
                    }

                }

            }

        }
    }


    sub export_sam_file {

        my ( $self, $sam_hash_ref ) = @_;
        my @header = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL NUMOFDIF ALIGNSCORE/;
        my @rows;

        foreach my $read ( sort {$a cmp $b} keys %{ $sam_hash_ref } ) {

            if ( $sam_hash_ref->{ $read }{ selected } eq "TRUE" ) {

                foreach my $pairs ( sort { $a cmp $b } keys %{ $sam_hash_ref->{ $read } } ) {

                    if ( $pairs =~ m/^\S+_[b|t]/ ) {

                        say join "\t", @{ $sam_hash_ref->{ $read }{ $pairs } }{ @header }  ;

                    }
                }
            }
        }

    }


    sub run {
        my ($self) = @_;
        my %bed6_hash;
        my %sam_hash;

        $self->bed6_to_hash( \%bed6_hash  );

        $self->clean_single_read( \%bed6_hash );

        %sam_hash = $self->sam_to_hash( $self->sam_file );

        $self->filter_reads( \%sam_hash, \%bed6_hash );

        $self->export_sam_file( \%sam_hash );

    }

    __PACKAGE__->meta->make_immutable;
}

use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

