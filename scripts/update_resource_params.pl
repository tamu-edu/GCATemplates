#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Temp qw/ tempfile tempdir /; # 'temp file';
use File::Copy 'move';

my $usage = qq{
Synopsis:
    This script will update the configurable Slurm resources in template scripts for time, ntasks-per-node, cpus-per-task and max-memory-per-node.
    Set the maximum CPU cores and memory to the maximum cores and memory availble on compute nodes.
        Only tools supporting multi-core processing will be updated with max cpus and memory values.
        Single core template scirpts will remain at 1 core.

Options:
    -n string       # cluster name
    -T int          # time in hours
    -K int          # tasks per node
    -C int          # cpus per task
    -M int          # max memory available on compute node in GB
    -m int          # min memory for single core jobs in GB
    -u string       # memory units to append to max memory or min memory. Default: G
    -h              # show usage

Example command:
    perl update_resource_params.pl -K 1 -C 48 -M 360

Configuration:
    A configuration line added at the bottom of each script is used to update tools that support multi-core processing
    and others that require the minimal amount of memory.
    The last line of the template script will need to contain the following format if you want to use update_resource_params.pl

    Example configs:
    #1:C:M	= set to 1 task: value provided by -C option for cores: and value provided by -M for memory
    #C:1:M	= set to -C tasks: 1 core: and -M for memory
    #1:1:m	= set to 1 task: 1 core: and -m for memory
};

my $cluster_name = 'hpc';
my $time_hours = 1;
my $ntasks_per_node = 1;
my $cpus_per_task = 12;
my $max_memory = 64;        # using units as set with $memory_units
my $min_memory = 2;         # using units as set with $memory_units
my $memory_units = 'G';

my $help;
GetOptions (
    # required
    "n|cluster_name=s" => \$cluster_name,   # default: hpc
    # optional
    "T|time_hours=i" => \$time_hours,           # default: 1
    "K|ntasks_per_node=i" => \$ntasks_per_node, # default: 1
    "C|cpus_per_task=i" => \$cpus_per_task,     # default: 12
    "M|max_memory=i" => \$max_memory,           # default in GB: 24
    "m|min_memory=i" => \$min_memory,           # used for single core jobs, default: $min_memory
    "u|memory_units=s" => \$memory_units,       # default: G
    "h|help" => \$help,
);

if (defined $help) {
    print "$usage\n";
    exit;
}

my @templates = glob("../templates/$cluster_name/run_*_${cluster_name}.sh");

foreach my $original_template_file (@templates) {
    my $old_file = get_content ($original_template_file);

    my ($tempfh, $tempfile) = tempfile();
    print $tempfh $old_file;
    close $tempfh;

    # Overwrite new file with old file
    move ($tempfile, $original_template_file) or die "Can't move \" $tempfile\"to \" $original_template_file\":$!";
}

sub get_content {
    my $file = shift;
    open my $fh, '<', $file or die "Can't open file \" $file\":$!";
    my $content = join('', <$fh>);

    # Slurm params that can be updated
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node

    if ($content =~ /\n#(\d+|C):(\d+|C):(\d+|M|m)$/) {
        my ($conf_ntasks, $conf_cpus, $conf_mem) = ($1,$2,$3) ;

        $content =~ s/--time=\d+:/--time=$time_hours:/;
        if ($conf_ntasks eq 'C') {
            $content =~ s/--ntasks-per-node=\d+/--ntasks-per-node=$cpus_per_task/;
        }
        elsif ($conf_ntasks eq 'K') {
            $content =~ s/--ntasks-per-node=\d+/--ntasks-per-node=$ntasks_per_node/;
        }
        $content =~ s/--cpus-per-task=\d+/--cpus-per-task=$cpus_per_task/ if $conf_cpus !~ /^\d+$/;
        $content =~ s/--mem=\d+\w/--mem=${max_memory}$memory_units/ if $conf_mem eq 'M';
        $content =~ s/--mem=\d+\w/--mem=${min_memory}$memory_units/ if $conf_mem eq 'm';
    }
    close $fh;
    return $content;
}

# example configuration parameters found on last line of each template file
# indicates recommended resources (ntasks, cpus-per-task, total-memory-per-node) for each templates script

#C:1:M
#1:C:M
#1:1:m
#1:2:m
#1:3:m
#2:1:m

# defined as:
# 1=keep at 1, do not update
# C=use $cpus_per_task
# M=use $max_memory
# m=use $min_memory
