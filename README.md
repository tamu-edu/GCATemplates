# GCATemplates
Genomic Computational Analysis Templates

This is a collection of HPC template job scripts designed to assist new HPC users in learning about submitting computational jobs from the command line. 

The templates contain example data sets and can be run immediately so users can see how the job submission process works.

After installing and configuring GCATemplates, run the command `gcatemplates` and a menu will be displayed allowing the user to select a category, task, tool and options for a template script which can then be copied to the current work directory.

## Installation:

`git clone https://github.com/tamu-edu/GCATemplates.git`

`cd GCATemplates/scripts`

`perl get_data.pl`

## Testing:

All scripts should be tested to verify they work on your cluster prior to making them available to users.

## Configuration:

The sample template job scripts in the repo are configured with 12 cores and 64GB of memory.
To update all multi-core scripts to use maximum cores available on your cluster, use the update_resource_params.pl tool.
For example for compute nodes with 48 cores and 360GB of available memory, do the following:

`cd GCATemplates/scripts`

`perl update_resource_params.pl -C 48 -M 360`

### Using your cluster names

The template script names have the cluster name at the end.

For example for a cluster named 'hpc', the template scripts will end with `_hpc.sh`

You can leave the names as _hpc.sh or change to your cluster name such as _gizmo.sh but you also have to edit the conf/templates.conf file with the new script names and bin/gcatemplates line that has the $system name variable.

### Adding new scripts

Put a copy of the new template script in the templates/hpc directory and add a line to the conf/templates.conf file.

### Adding scripts for a new cluster

In the templates directory, create a directory with the cluster name. For a cluster named Gizmo, create the directory `templates/gizmo` and save the template scripts in that directory and name the template scripts with suffix `_gizmo.sh`.
