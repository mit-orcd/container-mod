# container-mod

**container-mod** streamlines the process of pulling containers from public registries and automatically generating ready-to-use modulefiles. It is a versatile tool designed for use by HPC system administrators and group managers to create and manage modules accessible to all HPC users or group members. Additionally, regular users can leverage container-mod to create personal modulefiles for their individual workflows, enhancing efficiency and reproducibility in HPC environments.

## Overview

This script is designed to streamline the management of Singularity/Apptainer containers. It handles the following tasks:

1. Pulling container images: retrieve Singularity/Apptainer images from specified URIs.
2. Generating module files: create module files in Lmod format based on container images.
3. Creating Executable Wrappers: Generate wrapper bash scripts for each program or command provided by the target application. These lightweight scripts simplify access to containerized tools by mimicking the behavior of native executables, allowing users to run commands as if the application were installed locally. 
4. Pipeline Mode: Automatically pull images, generate module files, and create executables in one command.

The script also provides flexible options, including specifying custom output directories, forcing overwrites of existing files, and updating repository files with new version information. Additionally, for containers that include **Python** and **ipykernel**, it can also generate a Jupyter kernel file, enabling users to run the container within Jupyter Lab or Notebook on the Open OnDemand platform.

## Prerequisites

- Singularity or Apptainer: This script requires either Singularity or Apptainer to be installed and available in your environment.
- Lmod: For generating module files, the Lmod environment module system should be available.

## Usage

### Subcommands

- pull <URI>: Pulls a container image from the provided URI.
- module <URI>: Generates a module file for the container.
- exec <URI>: Creates wrapper bash scripts for the container’s programs.
- pipe <URI>: Executes the pipeline that pulls the image, generates a module file, and creates executables in one step.

### Options

 - -d, --dir DIR: Specify the output directory for images, module files, and executables. Defaults to the current directory.
- -f, --force: Force overwrite of existing module files, or executables. Default is to skip existing files.
-  -j, --jupyter: Generate Jupyter kernels for the specified URIs. With the jupyter kernel, users can run containers within Jupyter Lab or Jupyter Notebook. The prerequisite is that **ipykernel** must be installed within the container.
- -m, --moduledir DIR: Specify the directory that stores module files that can be used as template. Defaults to modulefiles.
- -u, --update: If set, the repository app file will be updated with new version information.
- -p, --personal: Create personal module files in users' `$HOME/privatemodules` (default is no).
- -h, --help: Display this help message and exit.

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

- images: singularity images will be stored in this folder.
- repos: The application information databases will be copied from the central copy to here. If users want to pull applicaitons that are missing in repos, they have to create one for these new applications.
- tools: bash wrapper for the supported executables will be saved here.

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
