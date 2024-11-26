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
    •	-p, --personal: Create personal module files in the privatemodules directory (default is no).
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

## Personal mode

If regular users want to pull containers and create modules for personal usage, users can run it in personal mode by adding `-p` or `--personal`.
This will create a folder named `container-apps` in users' `$HOME`. Within `container-apps`, there are three folders:

    - 1. images: singularity images will be stored in this folder.
    - 2. repos: The application information databases will be copied from the central copy to here. If users want to pull applicaitons that are missing in repos, they have to create one for these new applications.
    - 3. tools: bash wrapper for the supported executables will be saved here.

The modulefiles will be created in `$HOME/privatemodules`.

### Example:

```
container-mod pipe -p docker://staphb/bowtie2:2.5.4
```

Once completes, users can load personal modules and run the application:

```
$ module load use.own
$ module load bowtie2/2.5.4
$ bowtie2 --help
```


## Jupyter kernel

For containerized applications that are able to run on Jupyter notebook/lab, container-mod is also able to create a jupyter kernel. The kernel will be stored in users' $HOME under `$HOME/.local/share/jupyter/kernels`. To create jupyter kernels, users need to run the script with `-j` or `--jupyter`.
### Example
```
container-mod pipe -p -j docker://tensorflow/tensorflow:2.18.0-jupyter
```