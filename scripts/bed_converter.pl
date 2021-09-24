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

package MyApp::Bed6ToBed12 {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use namespace::autoclean;
    use Data::Printer;
    use List::MoreUtils qw( zip );

    command_long_description q[Converting bed6 reads file to a bed12 representing fragments];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
    );


    sub bed6_to_hash {

        my ( $self, $hash_ref ) = @_;
        my @header = qw/chr start end read_id score strand/;

        open( my $fh, '<', $self->input_file ) or die "Could not open file '$self->input_file' $!";
        while ( my $row = <$fh> ) {
            chomp $row;
            my @split_row = split( "\t", $row );

            if ( $split_row[ 4 ] >= 20 ) {

                my $read_id;
                my $read_type;
                if ( $split_row[ 3 ] =~ m/(\S+)_(\S+)\/.*/ ) {
                    $read_id   = $1;
                    $read_type = $2;
                }

                my %bed6_row_hash = zip @header, @split_row;

                $hash_ref->{ $read_id }{ $read_type } = \%bed6_row_hash

            }

        }

    }


    sub check_pairs {
        my ( $self, $hash_ref ) = @_;

        foreach my $pair ( keys %{ $hash_ref } ) {

            if ( !$hash_ref->{ $pair }{ 'bait' } or !$hash_ref->{ $pair }{ 'target' } ) {
                delete( $hash_ref->{ $pair } );
            }

        }

    }


    sub bed6_hash_to_bed12 {
        my ( $self, $hash_ref ) = @_;
        my %bed12_hash;


        foreach my $pair ( sort { $a cmp $b } keys %{ $hash_ref } ) {
            my $second_block_start;

            #check strand
            $bed12_hash{ $pair }{ 'strand' } = $hash_ref->{ $pair }{ 'bait' }{ 'strand' };

            if ( $bed12_hash{ $pair }{ 'strand' } eq "+" ) {
                $bed12_hash{ $pair }{ 'start' }      = $hash_ref->{ $pair }{ 'bait' }{ 'start' };
                $bed12_hash{ $pair }{ 'end' }        = $hash_ref->{ $pair }{ 'target' }{ 'end' };
                $bed12_hash{ $pair }{ 'thickStart' } = $hash_ref->{ $pair }{ 'bait' }{ 'start' };
                $bed12_hash{ $pair }{ 'thickEnd' }   = $hash_ref->{ $pair }{ 'target' }{ 'end' };
                $bed12_hash{ $pair }{ 'blockSizes' } = "1,10";

                $second_block_start = abs( ( ( $hash_ref->{ $pair }{ 'target' }{ 'end' } - 10 ) - $hash_ref->{ $pair }{ 'bait' }{ 'start' } ) );

            }
            else {
                $bed12_hash{ $pair }{ 'start' }      = $hash_ref->{ $pair }{ 'target' }{ 'start' };
                $bed12_hash{ $pair }{ 'end' }        = $hash_ref->{ $pair }{ 'bait' }{ 'end' };
                $bed12_hash{ $pair }{ 'thickStart' } = $hash_ref->{ $pair }{ 'target' }{ 'start' };
                $bed12_hash{ $pair }{ 'thickEnd' }   = $hash_ref->{ $pair }{ 'bait' }{ 'end' };
                $bed12_hash{ $pair }{ 'blockSizes' } = "10,1";

                $second_block_start = abs( ( $hash_ref->{ $pair }{ 'target' }{ 'start' } - ( $hash_ref->{ $pair }{ 'bait' }{ 'end' } - 1 ) ) );

            }

            $bed12_hash{ $pair }{ 'blockStarts' } = "0," . $second_block_start;
            $bed12_hash{ $pair }{ 'chrom' }       = $hash_ref->{ $pair }{ 'bait' }{ 'chr' };
            $bed12_hash{ $pair }{ 'name' }        = $pair;
            $bed12_hash{ $pair }{ 'score' }       = 1;
            $bed12_hash{ $pair }{ 'itemRgb' }     = 0;
            $bed12_hash{ $pair }{ 'blockCount' }  = 2;


        }

        return ( %bed12_hash );

    }


    sub export_bed12 {

        my ( $self, $hash_ref ) = @_;
        my @header_bed12 = qw/chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts /;

        foreach my $name ( keys %{ $hash_ref } ) {

            say join "\t", @{ $hash_ref->{ $name } }{@header_bed12}

        }
    }


    sub run {
        my ( $self ) = @_;
        my %bed_ref;

        $self->bed6_to_hash( \%bed_ref );

        $self->check_pairs( \%bed_ref );

        my %bed12 = $self->bed6_hash_to_bed12( \%bed_ref );

        $self->export_bed12(\%bed12)


    }

    __PACKAGE__->meta->make_immutable;
}

use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

