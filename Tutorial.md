# ANGSD-wrapper Tutorial
## Installation

Welcome! This is a short guide to population genetics analysis using ANGSD-wrapper. We will be using a test data set containing sequence from *Oryza sativa* and *Oryza glumaepetula*. First, we need to clone the ANGSD-wrapper repository. You need to have git installed to do this. Alternatively, you can download a zip file from the [releases page](https://github.com/arundurvasula/angsd-wrapper/releases) or use the download zip button on the home page of the repository.

### Dependencies

The basic dependencies for ANGSD-wrapper are `SAMTools`, `Wget`, and `Git`. Most Linux distributions have `Wget` and `Git` installed by default; users will need to download `SAMTools`, and it's dependency [`HTSlib`](https://github.com/samtools/htslib) from its [GitHub page](https://github.com/samtools/samtools).

### Note: Mac Users Have Special Installation Requirements

You will need to install all three of the basic dependencies to run, as well as the GNU Scientific Library. We recommend using [Homebrew](http://brew.sh/) to manage the installation process. If you use Homebrew, you can install all of the required dependencies from anywhere in your terminal with the following commands:

```shell
brew install git
brew install samtools
brew install wget
brew install gsl
```

### Downloading and Installing ANGSD-wrapper

We'll use `Git` to download ANGSD-wrapper. To do this, type the following commands:

```shell
git clone https://github.com/mojaveazure/angsd-wrapper.git
cd angsd-wrapper
```

ANGSD-wrapper comes with its own version of ANGSD to prevent compatibility breaking changes in ANGSD from affecting ANGSD-wrapper and comes with a few other programs. In order to compile these programs, you must run the setup routine.

```shell
./angsd-wrapper setup please
source ~/.bash_profile
```

This will download and install ANGSD, ngsAdmix, ngsTools, and ngsF. All of these programs are downloaded to the `dependencies` directory. In addition, it will also download and set up a directory with test data. These data are located in the `iplant` directory. Finally, ANGSD-wrapper will be installed system-wide so that it can be used from any working directory.

In the `iplant` directory, there are 12 BAM and BAI files (`[0-11].sub.bam` and `[0-11].sub.bam.bai`), serving as samples and their indeies, ancestral (`ancestral.merid_japonica_chr.fa`) and reference (`reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa`) sequences, a file with inbreeding coefficients (`InbreedingCoefficients.txt`), a list of regions (`regions.txt`) to be analyzed, and a file with sample names (`SampleNames.txt`).

ANGSD-wrapper has many different routines, or wrappers, that it can perform on a given dataset; we will be working using the Site Frequency Spectrum (SFS), Thetas Estimator, Admixture Analysis, and Principal Component Analysis (PCA) routines for this tutorial. We will also be graphing our results using ~~a Shiny web app~~ R code that will generate PDFs with plots for all outputs we generate. All analysis should be done using a supercomputer-like device, at least 32 GB of RAM~~, and all graphing should be done using a computer with a graphical user interaface~~. ~~If you have access to a supercomputer cluster, we recommend setting up ANGSD-wrapper on both the cluster for analysis and local machine for graphing.~~ To see all available wrappers, run `angsd-wrapper` without any arguments.

## Configuring ANGSD-wrapper with the `Common_Config` file

angsd-wrapper uses configuration files to figure out where the data is and what options should be passed to ANGSD and other dependencies. There is one configuration file per wrapper included with `angsd-wrapper`, as well as a common configuration file that can be used by multiple wrappers. All of these are located in the `Configuration_Files` directory; we recommend copying this directory to another directory so that there is always a clean copy of the configuration files availabl. In this case, we will copy the `Configuration_Files` directory into the `iplant` directory using the following command while being in the `angsd-wrapper` directory:

```shell
cp -r Configuration_Files/ iplant/
```

Now, let's go into the `iplant` directory and figure out the full path to this directory.

```shell
cd iplant/
pwd
```

This will output a string that starts with `/home/`, go ahead and copy everything following that second forward slash. For example, if we get `/home/software/angsd-wrapper/iplant` as our output, we only need `/software/angsd-wrapper/iplant`

Now, we'll go find our configuration files in the `Configuration_Files` directory:

```shell
cd Configuration_Files/
```

Because we're using multiple wrappers in this tutorial, we'll use the `Common_Config` file to hold variables that will be used across all methods. Open `Common_Config` in your favorite text editor, such as Vim or Emacs.

First, we need to define a list of samples. On line 10 of `Common_Config`, there's a place to define this sample list. If we remember back in our `iplant` directory, our sample list is called `SampleNames.txt`

So, to tell ANGSD-wrapper where our sample list is, using our exapmle of `iplant` being at `/home/software/angsd-wrapper/iplant`, we would make sure line 10 looks like this:

```shell
SAMPLE_LIST=${HOME}/software/angsd-wrapper/iplant/SampleNames.txt
```

### Note: BAM Files **MUST** Have an `@HD` Header Line

Some programs, when generating BAM files, will not include the `@HD` header line. To see if you have this line, use `SAMTools` to check the header for your BAM files:

```shell
samtools view -H <name of header> | head
```

The `@HD` header line should be the first line that pops up; if you don't see it, this [Gist](https://gist.github.com/mojaveazure/d194c4705642eecf8437) will add one for you.

Adjust the `/software/angsd-wrapper/iplant` part to whatever you copied from your output.

Next, we need our list of inbreeding coefficients. This was called `InbreedingCoefficients.txt`, we tell ANGSD-wrapper where this file is on line 13 of our `Common_Config` file:

```shell
SAMPLE_INBREEDING=${HOME}/software/angsd-wrapper/iplant/InbreedingCoefficients.txt
```

Lines 16 and 19 ask for our ancestral and reference sequences. These are `ancestral.merid_japonica_chr.fa` and `reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa`, respectively. In the `Common_Config` file, we'd enter the following on their respective lines:

```shell
ANC_SEQ=${HOME}/software/angsd-wrapper/iplant/ancestral.merid_japonica_chr.fa
REF_SEQ=${HOME}/software/angsd-wrapper/iplant/reference.Oryza_sativa.IRGSP-1.0.23.dna.genome_chr.fa
```

Now we need to set up our outdirectory structure. We use two variables to define this: `PROJECT` and `SCRATCH`. All output files will be placed in `$SCRATCH/$PROJECT/<name_of_program>`; for example, if we set `SCRATCH` to be "`${HOME}/scratch`" and `PROJECT` to be "Rice" and calculated a site frequency spectrum, our outdirectory would be `${HOME}/scratch/Rice/SFS`.

If we used the same `SCRATCH` and `PROJECT` assignments and estimated Thetas, our outdirectory would be `${HOME}/scratch/Rice/Thetas`. The outdirectory structure is generated automatically, making any directory within the structure that doesn't already exist, so there's need to make these directories before hand.

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

Now we're ready to run ANGSD-wrapper to calculate a site frequency spectrum, estimate Thetas, perform an admixture analysis, and run a principal component analysis. Close out of `Common_Config`, and be sure to save your changes. We're going to stay in the `Configuration_Files` directory for now.

## Site Frequency Spectrum

Each wrapper function has its own configuration file associated with it. To run the site frequency spectrum, we need the `Site_Frequency_Spectrum_Config` file. Open this up in your favorite text editor.

Each wrapper-specific configuration file is split into three parts: the `COMMON` definition, the 'not-using-common' section, and the wrapper-specifc variables section. If a wrapper utilizes the `Common` definition, it will always be on line 10. The 'not-using-common' section is blocked off by 

#Theta calculation

Running estimation of theta is very simple now that we've gotten everything set up. We need to modify one line in the `THETAS.sh` script:

    TAXON=test

If you would like to change the window size and step amount, you can modify the following variables:

    WIN=100
    STEP=50

You can also turn off sliding windows with:

    SLIDING_WINDOW=false

Great! Now we can run the script:

    sbatch -p bigmemm scripts/THETAS.sh scripts/THETAS.conf

We can visualize this the same way as before with the site frequency spectrum. The file we are interested in is called test_Diversity.thetas.gz.pestPG:

    $scp name@server.address.edu:~/rilab/aw-tutorial/results/test_Diversity.thetas.gz.pestPG ./

Once you upload the file into the web app, you should get two graphs that look like this:
![](http://i.imgur.com/gZemo59.png)
![](http://i.imgur.com/ZjfQhhZ.png)

You can view these graphs interactively by using the controls on the side. You can switch between estimators of theta and neutrality statistics with the drop down menus and zoom in on regions of interest with the Base start and end position number fields. You can also add lowess curves to the graphs and include gene annotations by uploading a GFF3 file that corresponds to your species.

## Admixture analysis

In order to carry out admixture analysis, we have to compile the ngsAdmix program.

    cd ngsPopGen
    g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

Then we just have to change one variable in the `ADMIX.conf` file.

    TAXON=test

The variable `K` specifies the ending point for the number of ancestral populations to run admixture with (starting ith 2). For example, if `K=5`, admixture will be run from K={2,3,4,5}. Now we can submit it like this:

    sbatch -p bigmemm scripts/ADMIX.sh scripts/ADMIX.conf

Once this has finished running, you can copy over the results (`*.qopt`) and view them in the graphing application.

    $scp name@server.edu:~/rilab/aw-tutorial/results/test.3.qopt ./

The results should look something like this:

![](http://i.imgur.com/hTYhnTo.png)

##PCA

Running a PCA is pretty straightforward. The only variable you need to change is the `TAXON` variable:

    TAXON=test

And we can run the analysis with

    sbatch -p bigmemm scripts/PCA.sh scripts/PCA.conf

We can visualize the result by copying over the `.covar` file:

    $scp name@server.edu:~/rilab/aw-tutorial/results/test_geno.covar ./

In the graphing application, switch to the PCA tab and upload the file. You should get something like this:

![](http://i.imgur.com/RMOheV1.png)
