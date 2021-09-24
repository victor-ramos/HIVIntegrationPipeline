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

package MyApp::GetReads {
    use feature qw(say);
    use MooseX::App::Command;
    extends 'MyApp';    # inherit log
    use MooseX::FileAttribute;
    use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
    use namespace::autoclean;
    use Data::Printer;
    use Capture::Tiny ':all';

    command_short_description q[This command is awesome];
    command_long_description q[This command is so awesome, yadda yadda yadda];

    has_directory 'fastq_dir' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );

    has_directory 'analysis_dir' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );


    sub initialize_dict {
        my ($self, $fastq_dict_ref) = @_;

        my @bait_fastq_files = glob( $self->analysis_dir->stringify . "/pre_selected_reads/*/*_bait_L*\.fq" );

        foreach my $fastq (@bait_fastq_files) {

            my $sample_name;
            if ($fastq =~ m/(\S+)pre_selected_reads\/(\S+)\/.*/) {
                $sample_name = $2;
            }

            $fastq_dict_ref->{$sample_name}{"raw_reads"} = 0;
            $fastq_dict_ref->{$sample_name}{"pre_selected_reads"} = 0;
            $fastq_dict_ref->{ $sample_name }{ "real_bait_target_reads" } = 0;
            $fastq_dict_ref->{$sample_name}{"reads_after_trimming"} = 0;
            $fastq_dict_ref->{ $sample_name }{ "human_reads" } = 0;
            $fastq_dict_ref->{ $sample_name }{ "hiv_reads" } = 0;
            $fastq_dict_ref->{ $sample_name }{ "filtered_human_reads" } = 0;

        }

    }


    sub get_raw_reads {

        my ($self, $fastq_dict_ref) = @_;

        my @fastq_files = glob( $self->fastq_dir->stringify . "/*.fastq.gz" );

        # Getting read count from raw fastq
        foreach my $fastq (@fastq_files) {

            if ( $fastq =~ m/(\S+)\/(\S+)_L.*/ ) {
                my $sample_name = $2;
                my $cmd = "gunzip -c ". $fastq . " | wc -l";

                my ( $stdout, $stderr, $exit ) = capture {
                    system( $cmd );
                };

                $fastq_dict_ref->{$sample_name}{"raw_reads"} = ( $stdout / 4 );
            }

        }

    }


    sub get_pre_selected_reads {

        my($self, $fastq_dict_ref) = @_;

        my @bait_fastq_files = glob( $self->analysis_dir->stringify . "/pre_selected_reads/*/*_bait_L*\.fq" );
        #my @target_fastq_files = glob( $self->analysis_dir->stringify . "/pre_selected_reads/*/*_target_L*" );
        #my @fastq_files = (@bait_fastq_files, @target_fastq_files);

        foreach my $fastq (@bait_fastq_files) {

            my $sample_name;
            if ($fastq =~ m/(\S+)pre_selected_reads\/(\S+)\/.*/) {
                $sample_name = $2;
            }

            my $cmd = "grep '\@M' ". $fastq . " | wc -l";

            my ( $stdout, $stderr, $exit );
            ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            $fastq_dict_ref->{$sample_name}{"pre_selected_reads"} = $stdout + 0;

        }

    }


    sub  get_real_baits_count {
        my ( $self, $fastq_dict_ref ) = @_;

        my @recovered_ltr_reads = glob( $self->analysis_dir->stringify . "/recover_ltr_reads/*/*.list" );
        foreach my $list ( @recovered_ltr_reads ) {

            my $sample_name;
            if ( $list =~ m/(\S+)recover_ltr_reads\/(\S+)\/.*/ ) {
                $sample_name = $2;
            }

            my $cmd =  "wc -l " . $list;

            my ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            my $read_number;
            if ( $stdout =~ m/(^\d+)\s.*/ ) {
                $read_number = $1;
            }

            $fastq_dict_ref->{ $sample_name }{ "real_bait_target_reads" } = ( $read_number );

        }

    }


    sub get_reads_after_trimming {

        my ($self, $fastq_dict_ref) = @_;

        my @bait_fastq_files = glob( $self->analysis_dir->stringify . "/bbduk_trim/*/*_bait_*" );
        foreach my $fastq (@bait_fastq_files) {

            my $sample_name;
            if ($fastq =~ m/(\S+)bbduk_trim\/(\S+)\/.*/) {
                $sample_name = $2;
            }

            my $cmd = "wc -l " . $fastq;

            my ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            if ( $stdout =~ m/(^\d+)\s\S+/ ) {
                $stdout = $1;
            }

            $fastq_dict_ref->{$sample_name}{"reads_after_trimming"} = ( $stdout / 4 );

        }


    }


    sub reads_aligned_human_genome {
        my ($self, $fastq_dict_ref) = @_;

        my @recovered_ltr_reads = glob( $self->analysis_dir->stringify . "/filtered_human_reads_alignment/*/*_proper_mapped_bait*.sam" );
        foreach my $list ( @recovered_ltr_reads ) {

            my $sample_name;
            if ( $list =~ m/(\S+)filtered_human_reads_alignment\/(\S+)\/.*/ ) {
                $sample_name = $2;
            }

            my $cmd =  "wc -l " . $list;


            my ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            my $read_number;
            if ( $stdout =~ m/(^\d+)\s.*/ ) {
                $read_number = $1;
            }

            $fastq_dict_ref->{ $sample_name }{ "human_reads" } = ( $read_number / 2);

        }

    }


    sub reads_aligned_hiv_genome {
        my ($self, $fastq_dict_ref) = @_;

        my @recovered_ltr_reads = glob( $self->analysis_dir->stringify . "/reads_alignment/*/*_hiv.bam" );
        foreach my $list ( @recovered_ltr_reads ) {

            my $sample_name;
            if ( $list =~ m/(\S+)reads_alignment\/(\S+)\/.*/ ) {
                $sample_name = $2;
            }

            my $cmd =  "samtools view -f 3 " . $list . " | wc -l";


            my ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            my $read_number;
            if ( $stdout =~ m/(^\d+)\s.*/ ) {
                $read_number = $1;
            }

            $fastq_dict_ref->{ $sample_name }{ "hiv_reads" } = ( $read_number / 2);

        }

    }


    sub filtered_reads {
        my ( $self, $fastq_dict_ref ) = @_;
        my @filtered_reads_bed;
        my $pattern = "filtered";
        my @filtered_reads = glob( $self->analysis_dir->stringify . "/human_bed_and_hostspot/*/*bait_target.bed" );
        @filtered_reads_bed = grep { !/$pattern/ } @filtered_reads;

        foreach my $bed ( @filtered_reads_bed ) {

            my $sample_name;
            if ( $bed =~ m/(\S+)human_bed_and_hostspot\/(\S+)\/.*/ ) {
                $sample_name = $2;
            }

            my $cmd =  "wc -l " . $bed;

            my ( $stdout, $stderr, $exit ) = capture {
                system( $cmd );
            };

            my $read_number;
            if ( $stdout =~ m/(^\d+)\s.*/ ) {
                $read_number = $1;
            }

            $fastq_dict_ref->{ $sample_name }{ "filtered_human_reads" } = ( $read_number );

        }

    }


    sub export_table {
        my ($self, $fastq_dict_ref) = @_;

        my @header = ('raw_reads', 'pre_selected_reads', 'real_bait_target_reads', 'reads_after_trimming','human_reads', 'hiv_reads', 'filtered_human_reads');


        say join "\t", ("sample", @header);
        foreach my $sample ( keys %{$fastq_dict_ref} ) {

            my @rows;
            foreach my $column ( @header ) {

                push @rows, $fastq_dict_ref->{$sample}{$column};

            }


            say join "\t", ($sample, @rows);
        }


    }


    sub run {
        my ($self) = @_;

        my %fastq_dict;

        $self->initialize_dict( \%fastq_dict );

        $self->get_raw_reads( \%fastq_dict );

        $self->get_pre_selected_reads( \%fastq_dict );

        $self->get_real_baits_count( \%fastq_dict  );

        $self->get_reads_after_trimming( \%fastq_dict  );

        $self->reads_aligned_human_genome( \%fastq_dict  );

        $self->reads_aligned_hiv_genome( \%fastq_dict  );

        $self->filtered_reads( \%fastq_dict  );

        $self->export_table( \%fastq_dict  );

    }

    __PACKAGE__->meta->make_immutable;
}

use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly
MyApp->new_with_command->run();

