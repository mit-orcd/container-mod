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

- `pull <URI>`: Pulls a container image from the provided URI.
- `module <URI>`: Generates a module file for the container.
- `exec <URI>`: Creates wrapper bash scripts for the container’s programs.
- `pipe <URI>`: Executes the pipeline that pulls the image, generates a module file, and creates executables in one step.

### Options

 - `-d, --dir DIR`: Specify the output directory for images, module files, and executables. Defaults to the current directory.
- `-f, --force`: Force overwrite of existing module files, or executables. Default is to skip existing files.
- `-j, --jupyter`: Generate Jupyter kernels for the specified URIs. With the jupyter kernel, users can run containers within Jupyter Lab or Jupyter Notebook. The prerequisite is that **ipykernel** must be installed within the container.
- `-m, --moduledir DIR`: Specify the directory that stores module files that can be used as template. Defaults to modulefiles.
- `-u, --update`: If set, the repository app file will be updated with new version information.
- `-p, --personal`: Create personal module files in users' `$HOME/privatemodules` (default is no).
- `-h, --help`: Display this help message and exit.

## Repository database
In the repos folder, each scientific application is represented by an individual information file, which serves as the basis for generating the corresponding module file. Each info file requires three fields: **Description**, **Home Page**, and **Programs**. The **Description** and **Home Page** fields provide metadata for the module file, while **Programs** is a comma-delimited list of commands or executables provided by the application. The workflow uses this list to generate bash wrappers for each command. Additionally, an optional **version** field helps track the URI of the pulled containers for reproducibility.

```
Description: Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.
Home Page: https://github.com/BenLangmead/bowtie2
Programs: bowtie2,bowtie2-build,bowtie2-inspect

version("2.5.4", uri="docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4")
version("2.5.1", uri="docker://quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0")
version("2.4.2", uri="docker://quay.io/biocontainers/bowtie2:2.4.2--py38hc2f83ea_2")
```

## Examples

### Pull an image

```
$ container-mod pull docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Generate a lmod module file

```
$ container-mod module docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Create bash wrappers for executables

```
$ container-mod exec docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Run the whole pipeline (pull, module and exec)

```
$ container-mod pipe docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

### Force overwrite existing files

If output files already exist in the specified output folder, the script will skip generating new files by default. To overwrite existing files and regenerate the output, users can include the `-f` or `--force` option.

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
