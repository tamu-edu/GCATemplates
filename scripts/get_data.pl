#!/usr/bin/env perl

use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /; # 'temp file';
use File::Copy 'move';
use File::Spec;

use Getopt::Long qw(:config no_ignore_case);

my $usage = qq{
Synopsis:
    This script will download data needed for GCATemplates template scripts to run directly.
    Run this script from the scripts directory which usually takes less than 5 minutes to complete.
        cd GCATemplates/scripts
        perl get_data.pl
};

my $help;
my $data_dir = "../data";
my $cluster_name = 'hpc';

GetOptions (
    # optional
    "n|cluster_name=s" => \$cluster_name,   # default: hpc
    "d|data_dir=s" => \$data_dir,           # default: GCATemplates/data/ inside install directory
    "h|help" => \$help,
);

if (defined $help) {
    print "$usage\n";
    exit;
}

# genomic E. coli
my $genomic_e_coli = qx{
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR237/097/SRR23746297/SRR23746297_1.fastq.gz -o $data_dir/e_coli/seqs/SRR23746297_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR237/097/SRR23746297/SRR23746297_2.fastq.gz -o $data_dir/e_coli/seqs/SRR23746297_2.fastq.gz
};

# transcriptomic E. coli
my $ecoli_data_1 = qx{
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/070/ERR10828770/ERR10828770_1.fastq.gz -o $data_dir/e_coli/seqs/ERR10828770_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/070/ERR10828770/ERR10828770_2.fastq.gz -o $data_dir/e_coli/seqs/ERR10828770_2.fastq.gz
};

# metagenomic
my $metagenomic_data = qx{
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR15366205/SRR15366205_1.fastq.gz -o $data_dir/metagenomic/SRR15366205_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR15366205/SRR15366205_2.fastq.gz -o $data_dir/metagenomic/SRR15366205_2.fastq.gz
};

# for macs2 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3868217/
my $chip_seq = qx{
curl -L http://cistrome.dfci.harvard.edu/MACSNatureProtocol/UW_K562_H3K4me3.tar.gz -o $data_dir/chipseq/UW_K562_H3K4me3.tar.gz
cd $data_dir/chipseq && tar xzf UW_K562_H3K4me3.tar.gz && rm UW_K562_H3K4me3.tar.gz && cd -
};

my @templates = glob("../templates/$cluster_name/run_*_${cluster_name}.sh");
my $data_dir_abs_path = File::Spec->rel2abs( $data_dir );
$data_dir_abs_path =~ s,scripts/../,,;

foreach my $original_template_file (@templates) {
    my $old_file = get_content ($original_template_file);

    my ($tempfh, $tempfile) = tempfile();
    print $tempfh $old_file;
    close $tempfh;

    move ($tempfile, $original_template_file) or die "Can't move \" $tempfile\"to \" $original_template_file\":$!";
}

sub get_content {
    my $file = shift;
    open my $fh, '<', $file or die "Can't open file \" $file\":$!";
    my $content = join('', <$fh>);

    $content =~ s,DATA_DIR,$data_dir_abs_path,g;

    close $fh;
    return $content;
}
