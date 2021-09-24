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

package MyApp::FilterHumanAlignment {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use namespace::autoclean;
    use Data::Printer;

    command_short_description q[This command is awesome];
    command_long_description q[This command is so awesome, yadda yadda yadda];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );


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


    sub select_pairs {

        my ( $self, $sam_hash_ref ) = @_;

        foreach my $pair ( keys %{ $sam_hash_ref } ) {
            my $matches_count;
            my $mismatches_count;
            my $soft_clipping_count;

            foreach my $read ( keys %{ $sam_hash_ref->{ $pair } } ) {
                if ( $read =~ m/^\S+_bait$/ ) {
                    my $cigar = $sam_hash_ref->{ $pair }{ $read }{ CIGAR };

                    if ( $cigar =~ m/^(\d+M)\S*/ ) {
                        $matches_count = $1;
                        $matches_count =~ s/M//g;

                        if ( $matches_count >= 30 ) {
                            $sam_hash_ref->{ $pair }{ selected } = "TRUE";
                        }
                        next;
                    }
                    elsif ( $cigar =~ m/^(\d+S)(\d+M)\S*/ ) {
                        $soft_clipping_count = $1;
                        $matches_count       = $2;

                        $soft_clipping_count =~ s/S//g;
                        $matches_count       =~ s/M//g;

                        if ( $soft_clipping_count >= 1 and $soft_clipping_count <= 3 and $matches_count >= 27 ) {
                            $sam_hash_ref->{ $pair }{ selected } = "TRUE";
                        }
                        next;
                    }
                    elsif ( $cigar =~ m/^(1)M([12])X(\d+)M\S*/  ) {
                        $matches_count       = $1 + $3;
                        $mismatches_count = $2;

                        if (  ( $matches_count - $mismatches_count )  >= 28 ) {
                            $sam_hash_ref->{ $pair }{ selected } = "TRUE";
                        }
                        next;
                    }
                }
            }
        }
    }


    sub export_sam_file {

        my ( $self, $sam_hash_ref ) = @_;
        my @header = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL NUMOFDIF ALIGNSCORE/;
        my @rows;

        foreach my $read ( keys %{ $sam_hash_ref } ) {

            if ( $sam_hash_ref->{ $read }{ selected } eq "TRUE" ) {

                foreach my $pairs ( sort { $a cmp $b } keys %{ $sam_hash_ref->{ $read } } ) {

                    if ( $pairs =~ m/^\S+_[b|t]/ ) {
                        say join "\t", @{$sam_hash_ref->{ $read }{ $pairs }}{@header};
                    }
                }
            }
        }

    }


    sub run {
        my ( $self ) = @_;

        my %sam_hash = $self->sam_to_hash( $self->input_file->stringify );

        $self->select_pairs( \%sam_hash );

        $self->export_sam_file( \%sam_hash );

    }

    __PACKAGE__->meta->make_immutable;
}

use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

