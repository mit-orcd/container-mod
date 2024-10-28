# container-mod

Workflow to pull containers from public registries and generate modulefiles
README for Container Management Script

## Overview

This script is designed to streamline the management of Singularity/Apptainer containers. It handles the following tasks:

    1.	Pulling Container Images: Retrieve Singularity/Apptainer images from specified URIs.
    2.	Generating Module Files: Create module files in Lmod format based on container images.
    3.	Creating Executables: Generate wrapper bash scripts that link to container programs.
    4.	Pipeline Mode: Automatically pull images, generate module files, and create executables in one command.

The script also supports options to specify output directories, force overwriting existing files, and update repository files with new version information.

## Prerequisites

    •	Singularity or Apptainer: This script requires either Singularity or Apptainer to be installed and available in your environment.
    •	Lmod: For generating module files, the Lmod environment module system should be available.

## Usage

### Subcommands

    •	pull <URI>: Pulls a container image from the provided URI.
    •	module <URI>: Generates a module file for the container.
    •	exec <URI>: Creates a wrapper bash script for the container’s programs.
    •	pipe <URI>: Executes a pipeline that pulls the image, generates a module file, and creates the executable in one step.

### Options

    •	-d, --dir DIR: Specify the output directory for images, module files, and executables. Defaults to the current directory.
    •	-f, --force: Force overwrite of existing module files, or executables. Default is to skip existing files.
    •	-m, --moduledir DIR: Specify the directory that stores module files that can be used as template. Defaults to modulefiles.
    •	-u, --update: If set, the repository app file will be updated with new version information.
    •	-h, --help: Display this help message and exit.

## Examples

### Pull an image:

```
$ container-mod pull docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Create a lmod module file

```
$ container-mod module docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Creating executables

```
$ container-mod exec docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Run the whole pipeline (pull, module and exec)

```
$ container-mod pipe docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Force overwrite existing files

```
$ container-mod pipe -f docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Update the repo database

```
$ container-mod pipe -u docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```
