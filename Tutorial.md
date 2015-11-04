# ANGSD-wrapper Tutorial
## Installation

Welcome! This is a short guide to population genetics analysis using ANGSD-wrapper. We will be using a test data set containing sequence from *Oryza sativa* and *Oryza glumaepetula*. First, we need to clone the ANGSD-wrapper repository. You need to have git installed to do this. Alternatively, you can download a zip file from the [releases page](https://github.com/arundurvasula/angsd-wrapper/releases) or use the download zip button on the home page of the repository.

### Dependencies

The basic dependencies for ANGSD-wrapper are `SAMTools`, `Wget`, and `Git`. Most Linux distributions have `Wget` and `Git` installed by default, however some will not; users will need to download `SAMTools`, and it's dependency [`HTSlib`](https://github.com/samtools/htslib) from its [GitHub page](https://github.com/samtools/samtools) or your package manager.

---
>### Note: Mac Users Have Special Installation Requirements

>You will need to install all three of the basic dependencies to run, as well as the GNU Scientific Library. We recommend using [Homebrew](http://brew.sh/) to manage the installation process. If you use Homebrew, you can install all of the required dependencies from anywhere in your terminal with the following commands:

>```shell
brew install git
brew install samtools
brew install wget
brew install gsl
```
---

### Downloading and Installing ANGSD-wrapper

We'll use `Git` to download ANGSD-wrapper. To do this, type the following commands:

```shell
git clone https://github.com/mojaveazure/angsd-wrapper.git
cd angsd-wrapper
```

ANGSD-wrapper comes with its own version of ANGSD to prevent compatibility breaking changes in ANGSD from affecting ANGSD-wrapper and comes with a few other programs. In order to compile these programs, you must run the following setup routine:

```shell
./angsd-wrapper setup please
source ~/.bash_profile
```

This will download and install ANGSD, ngsAdmix, ngsTools, and ngsF. All of these programs are downloaded to the `dependencies` directory. In addition, it will also download and set up a directory with test data. These data are located in the `iplant` directory. Finally, ANGSD-wrapper will be installed system-wide so that it can be used from any working directory. To make sure ANGSD-wrapper installed correctly, run `angsd-wrapper`, without the `./` that we used before.

In the `iplant` directory, there are 12 BAM and BAI files (`[0-11].sub.bam` and `[0-11].sub.bam.bai`), serving as samples and their indices, ancestral (`ancestral.merid_japonica_chr.fa`) and reference (`reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa`) sequences, a file with inbreeding coefficients (`InbreedingCoefficients.txt`), a list of regions (`regions.txt`) to be analyzed, and a file with sample names (`SampleNames.txt`).

ANGSD-wrapper has many different routines, or wrappers, that it can perform on a given dataset; we will be working with the Site Frequency Spectrum (SFS), Thetas Estimator, Admixture Analysis, and Principal Component Analysis (PCA) routines for this tutorial. To see all available wrappers, run `angsd-wrapper` without any arguments.

We will also be graphing our results using a Shiny web app. All analysis should be done using a supercomputer-like device, at least 32 GB of RAM, and all graphing should be done using a computer with a graphical user interaface. If you have access to a supercomputer cluster, we recommend setting up ANGSD-wrapper on both the cluster for analysis and local machine for graphing.

## Configuring ANGSD-wrapper with the `Common_Config` file

ANGSD-wrapper uses configuration files to figure out where the data is and what options should be passed to ANGSD and other dependencies. There is one configuration file per wrapper included with `angsd-wrapper`, as well as a common configuration file (`Common_Config`) that can be used by multiple wrappers. All of these are located in the `Configuration_Files` directory; we recommend copying this directory to another directory so that there is always a clean copy of the configuration files available. In this case, we will start in the `angsd-wrapper` directory; then we will copy the `Configuration_Files` directory into the `iplant` directory using the following command:

```shell
cp -r Configuration_Files/ iplant/
```

Now, let's go into the `iplant` directory and figure out the full path to this directory using `pwd`.

```shell
cd iplant/
pwd
```

This will output a string that starts with `/home/`; go ahead and copy everything following that second forward slash. For example, if we get `/home/software/angsd-wrapper/iplant` as our output, we only need `/software/angsd-wrapper/iplant`.

Now, we'll go find our configuration files in the `Configuration_Files` directory:

```shell
cd Configuration_Files/
```

Because we're using multiple wrappers in this tutorial, we'll use the `Common_Config` file to hold variables that will be used across all methods. Open `Common_Config` in your favorite [text editor](http://www.yolinux.com/TUTORIALS/LinuxTextEditors.html), such as [Vim](http://vim.wikia.com/wiki/Tutorial) or [Emacs](https://www.gnu.org/software/emacs/manual/).

First, we need to define a list of samples. On line 10 of `Common_Config`, there's a place to define this sample list. If we remember back in our `iplant` directory, our sample list is called `SampleNames.txt`

So, to tell ANGSD-wrapper where our sample list is, we will use our `iplant` example with the directory location being  `/home/software/angsd-wrapper/iplant`. We'll make sure line 10 looks like this:

```shell
SAMPLE_LIST=${HOME}/software/angsd-wrapper/iplant/SampleNames.txt
```

We use `${HOME}` to help ANGSD-wrapper find your files, if we used `/home`, some systems will error out, saying that the directory does not exist.

---
>### Note: BAM Files **MUST** Have an `@HD` Header Line

>Some programs, when generating BAM files, will not include the `@HD` header line. To see if you have this line, use `SAMTools` to check the header for your BAM files:

>```shell
>samtools view -H <name of BAM file> | head -1
>```

>The `@HD` header line should be the first line that pops up; if you don't see it, this [Gist](https://gist.github.com/mojaveazure/d194c4705642eecf8437) will add one for you.

---

Adjust the `/software/angsd-wrapper/iplant` part to whatever you copied from your output.

Next, we need our list of inbreeding coefficients. This is called `InbreedingCoefficients.txt`, to run this we tell ANGSD-wrapper where this file is on line 13 of our `Common_Config` file:

```shell
SAMPLE_INBREEDING=${HOME}/software/angsd-wrapper/iplant/InbreedingCoefficients.txt
```

Lines 16 and 19 ask for our ancestral and reference sequences. These are `ancestral.merid_japonica_chr.fa` and `reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa`, respectively. In the `Common_Config` file, we'd enter the following on their respective lines:

```shell
ANC_SEQ=${HOME}/software/angsd-wrapper/iplant/ancestral.merid_japonica_chr.fa
REF_SEQ=${HOME}/software/angsd-wrapper/iplant/reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa
```

Now we need to set up our outdirectory structure. We use two variables to define this: `PROJECT` and `SCRATCH`. All output files will be placed in `$SCRATCH/$PROJECT/<name_of_program>`; for example, if we set `SCRATCH` to be "`${HOME}/scratch`" and `PROJECT` to be "Rice" and calculated a site frequency spectrum, our outdirectory would be `${HOME}/scratch/Rice/SFS`.

If we used the same `SCRATCH` and `PROJECT` assignments and estimated Thetas, our outdirectory would be `${HOME}/scratch/Rice/Thetas`. The outdirectory structure is generated automatically, making any directory within the structure that doesn't already exist, so it is necessary to make these directories before hand.

Let's set `SCRATCH` to be "`${HOME}/scratch`" and `PROJECT` to be "Rice"; we define these two variables on lines 22, for `PROJECT`, and 27, for `SCRATCH`:

```shell
PROJECT=Rice
SCRATCH=${HOME}/Rice
```

Finally, we need to specify a regions file for ANGSD-wrapper. While we can run ANGSD-wrapper without a regions file, it becomes very computationally expensive and takes much longer. If you would like to generate a regions file, taking a random sample of all possible regions, this [Gist](https://gist.github.com/mojaveazure/d115bb25eeff3b2df9f9) will create a valid regions file for you.

We have a regions file in our `iplant` directory called `regions.txt`, let's tell ANGSD-wrapper where this is on line 31 of `Common_Config`

```shell
REGIONS=${HOME}/software/angsd-wrapper/iplant/regions.txt
```

Now we're ready to run ANGSD-wrapper to calculate a site frequency spectrum, estimate Thetas, perform an analyze admixture, and run a principal component analysis. Close out of `Common_Config`, and be sure to save your changes. We're going to stay in the `Configuration_Files` directory for now; we're going to need the full path for this directory, which we obtain with the following command:

```shell
pwd
```

Again, we'll get an output starting with `/home/` and we only need the part after the second forward slash. Using our directory structure from before, our output would be `/home/software/angsd-wrapper/iplant/Configuration_Files` and we need `/software/angsd-wrapper/iplant/Configuration_Files`

## Site Frequency Spectrum

Each wrapper function has its own configuration file associated with it. To run the site frequency spectrum, we need the `Site_Frequency_Spectrum_Config` file. Open this up in your favorite text editor.

---
>### Note: A Word About Configuration Files

>Each wrapper-specific configuration file is split into three parts: the `COMMON` definition, the 'not-using-common' section, and the wrapper-specifc variables section. If a wrapper utilizes the `Common` definition, it will always be on line 10. The 'not-using-common' section is blocked off by 94 hash marks (`#`). If you are not using the `Common_Config` file, please fill out the variable definitions in this section. Since we're using `Common_Config`, we can skip these lines. finally, the wrapper-specific section includes any other variable definitons as well as parameters for the specifc wrapper.

---

In `Site_Frequency_Spectrum_Config`, we need to tell ANGSD-wrapper where our `Common_Config` file is. This definition is on line 10:

```shell
COMMON=${HOME}/software/angsd-wrapper/iplant/Configuration_Files/Common_Config
```

Remember to adjust the `/software/angsd-wrapper/iplant/Configuration_Files` part for your own directory structure.

Most of the other variables and parameters are set up to run smoothly. For now, we're going to set `OVERRIDE` to be `true`, in case we run this another time and want updated results; we change this on line 53:

```shell
OVERRIDE=true
```

### Note: How ANGSD-wrapper Knows What to do

ANGSD-wrapper has several functions, or wrappers, built into it. These are predefined and very specific. The syntax for running ANGSD-wrapper is as follows:

```shell
angsd-wrapper <wrapper> <configuration file>
```

Where `<wrapper>` is one of the wrappers that ANGSD-wrapper can perform and `<configuration file>` is the full path to the configuration file we set up for it. To see a full list of wrappers that ANGSD-wrapper has and how to call them, run the following command:

```shell
angsd-wrapper
```
This will display a usage message with the wrappers ANGSD-wrapper has and how to call them. Capitalization and spelling are very important with ANGSD-wrapper; you **must** type out what you see in the usage message to get ANGSD-wrapper to run. Also, you don't have to use our presupplied configuration files with ANGSD-wrapper, but you *do* need to have all of the variable definitions that we have supplied.

Now, lets calculate a site frequency spectrum using ANGSD-wrapper:

```shell
angsd-wrapper SFS ./Site_Frequency_Spectrum_Config
```

Once this finishes, our output files will be in the outdirectory we specified, `${HOME}/scratch/Rice/SFS`, let's go there and look at our files:

```shell
cd ${HOME}/scratch/Rice/SFS/
ls
```

Here we see files. I don't know what they're called, but hopefully we figure this shit out.

We'll need the `_DerivedSFS` file for our Thetas estimation

## Thetas Estimation

Now, we need to go back to our `Configuration_Files` directory so we can set up ANGSD-wrapper to estimates Thetas values for us. We use the `cd` command to do this:

```shell
cd ${HOME}/software/angsd-wrapper/iplant/Configuration_Files/
```

Open up `Thetas_Config` in your favorite text editor. We have three variables we need to define in this configuration file. First, we need to tell ANGSD-wrapper where our `Common_Config` file is; this will be the same as what we put in our `Site_Frequency_Spectrum_Config` file. On line 10, we'll put:

```shell
COMMON=${HOME}/software/angsd-wrapper/iplant/Configuration_Files/Common_Config
```

Next, we need to specify our pest file. This file comes from our site frequency spectrum; in this case, our file is `_DerivedSFS`. We need to specify this on line 41 of `Thetas_Config`:

```shell
PEST=${HOME}/scratch/Rice/SFS/_DerivedSFS
```

Finally, let's set `OVERRIDE` to be `true` again; we do this on line 56:

```shell
OVERRIDE=true
```

Now, we can estimate Thetas values using ANGSD-wrapper; we do this with the following command:

```shell
angsd-wrapper Thetas ./Thetas_Config
```

Our output files will be in `${HOME}/scratch/Rice/Thetas`, let's go there and look at our files

```shell
cd ${HOME}/scratch/Rice/Thetas/
ls
```

Here we see more files, which I don't know what they are. They're super cool, though, and that one is my favorite. No, not that one, the one next to it. There you go. Isn't it awesome?

## Admixture Analysis

Let's go back to out `Configuration_Files` directory to set up our admixture analysis:

```shell
cd ${HOME}/software/angsd-wrapper/iplant/Configuration_Files/
```

We need to edit variables in `Admixture_Config` to tell ANGSD-wrapper where everything is for the admixture analysis. Open up `Admixture_Analysis` with your favorite text editor. On line 10, we need to specify where our `Common_Config` file is:

```shell
COMMON=${HOME}/software/angsd-wrapper/iplant/Configuration_Files/Common_Config
```

The only other variable we need to specify is our likelihood file. This comes from the site frequency spectrum and is called `*.beagle.gz`; remember, this is located in `${HOME}/scratch/Rice/SFS`. On line 26 of `Admixture_Config`, we need to tell ANGSD-wrapper where this likelihood file is:

```shell
LIKELIHOOD=${HOME}/scratch/Rice/SFS/.beagle.gz
```

Now, let's run the admixture analysis. This is done with the following command:

```shell
angsd-wrapper Admixture ./Admixture_Config
```

Our output files will be in `${HOME}/scratch/Rice/ngsAdmix`, let's go there an look at our files.

```shell
cd ${HOME}/scratch/Rice/ngsAdmix/
ls
```

Here, we see some more fantastic, unknown files. Let us ponder their existence and relish in their mystery.

## Principal Component Analysis

Let's go back to out `Configuration_Files` directory to set up our principal component  analysis (PCA):

```shell
cd ${HOME}/software/angsd-wrapper/iplant/Configuration_Files/
```

We only need to tell ANGSD-wrapper where our `Common_Config` file is, everythin else for the PCA is taken care of. We do this on line 10 of `Principal_Component_Analysis_Config`:

```shell
COMMON=${HOME}/software/angsd-wrapper/iplant/Configuration_Files/Common_Config
```

Now, let's run the PCA. We use the following command to do this:

```shell
angsd-wrapper PCA ./Principal_Component_Analysis_Config
```

Our output files will be in `${HOME}/scratch/Rice/PCA`, let's go there an look at our files.

```shell
cd ${HOME}/scratch/Rice/PCA/
ls
```

More files :D

## Graphing

ANGSD-wrapper comes with a visualization package, based off of Rstudio's Shiny platform. To use this, we need to be on a machine that has a graphical user interface (GUI) and a web browser. If you used ANGSD-wrapper on a high performance computing system, please transfer your files to another machine with a GUI so we can run utilize the visualization package. You may need to setup ANGSD-wrapper again.

The files we need for graphing are:
 - `_Derived_SFS`
 - `.pestPG`
 - All `.qopt` files
 - `_PCA.covar`

### Shiny Graphing

> **Hurry up _Chaochih_**
