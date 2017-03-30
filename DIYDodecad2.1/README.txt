Do-It-Yourself Dodecad 2.0
(DIYDodecad 2.0)

Copyright (c) 2011 Dienekes Pontikos

Contents

1) TERMS OF USE
2) WHAT'S NEW in 2.1
3) REQUIREMENTS
4) INSTRUCTIONS
5) TECHNICAL DETAILS
6) RESULT FILES
7) VISUALIZATION
8) FEEDBACK

1) TERMS OF USE

You are free to use this software and all the files contained within this 
directory, for any non-commercial reason, provided that you attribute them to 
Dienekes Pontikos and/or the Dodecad Ancestry Project, and that you provide a 
link to either:

http://dienekes.blogspot.com
http://dodecad.blogspot.com

The user accepts full responsibility for all consequences deriving from the use 
of this software.

2) WHAT'S NEW in 2.1

DIYDodecad v 2.1 allows incomplete genotype files to be used, i.e., genotype
files that do not include all expected SNP markers used in a calculator. This 
is useful to individuals having older genotype files from their testing 
companies, and allows the tool to be used with any type of genotype data. 

There is a minimum requirement of at least 100 usable SNPs, i.e., SNPs that are
in the genotype file and do not have no-calls.

DIYDodecad v 2.0 has several changes and new features compared to 1.0:

- Convergence is assessed by measuring the difference in admixture proportions 
in successive iterations; if any of them is greater than the "goal", iterations 
continue
- There are no longer any limitations on the number of SNPs/ancestral 
populations for those interested in developing their own calculators
- Results are written to files and not just to the screen

Most importantly:

- It is now possible to get admixture proportions for all chromosomes, all 
regions of a certain length, or a particular region of interest

3) REQUIREMENTS

DIYDodecad requires:
1. A Windows or 32bit/64bit Linux machine
2. The R software, which can be downloaded and installed from 
http://www.r-project.org/
3. Raw autosomal genotype data as can be obtained from 23andMe 
(http://www.23andme.com) or Family Finder Illumina 
(http://www.familytreedna.com/landing/family-finder.aspx)

CPU and memory requirements are modest, and results will, in most cases, be 
computed within a few minutes for the simple admixture analysis, and a few
hours at most for the more complicated segment analysis.

4) INSTRUCTIONS

The following instructions show how you can use DIYDodecad for the first time. 
You should have no problems if you follow them closely.

1. Unpack the contents of the DIYDodecad2.0.rar file to a directory in your 
computer. Henceforth, we will call this the "working directory".

If you don't know how to do this:

In Windows:

Double click on the DIYDodecad2.0.rar file to see if you already have software 
installed on your computer that can open it. Otherwise, you can use the free 
utility
7-zip (http://www.7-zip.org/) or WinRar (commercial, but functional even if you 
don't buy it http://www.rarlab.com/download.htm)

In Linux:

Enter the command:

unrar e DIYDodecad2.0.rar

If unrar is not installed, you need to install it; this will depend on your 
particular Linux installation; in Ubuntu, you can install it by:

sudo apt-get install unrar

2. Download and install R (http://www.r-project.org/) if you haven't done so 
already.

In Windows:

Select a mirror that is close to you:

http://cran.r-project.org/mirrors.html

and then click on "base" to get to the installation file.

In Linux:

Again, this will depend on your particular Linux installation; in Ubuntu you 
can install it by:

sudo apt-get install r-base-core

Alternatively, you can also go to a mirror close to you:

http://cran.r-project.org/mirrors.html

and then follow the instructions to install it for your operating system.

3. Download and unzip your genotype data file from 23andMe or FamilyTreeDNA 
into the working directory. This should be a file ending in .txt for 23andMe 
and .csv for
FamilyTreeDNA. If you don't know how to uncompress your downloaded data, you 
can use the same software as in step #1 above.

From now on, I will suppose that your genotype file is named 'johndoe.txt' or 
'johndoe.csv'

4. Launch R. In Linux you just type:

R

In Windows, R will be listed in your Program Files or All Programs from the 
Start Menu, and you can start it like any other Windows Program.

Once R is running, it will give you a command prompt where you can enter 
commands. First, you must change the directory to your working directory.

You can do this from the File -> Change dir menu in Windows, or by using the 
setwd command in either Windows or Linux, e.g.:

setwd('/home/johndoe/dodecad')
setwd('c:\\users\\johndoe\\Dodecad')

Then, enter:

source('standardize.r')

This loads a small program that will convert your data from the 
company-specific format to a common format in the next step.

5. At the R prompt, enter: 

a. If you have 23andMe data (either v2 or v3 chip): 

standardize('johndoe.txt', company='23andMe')

b. if you have Family Finder data (Illumina chip only):

standardize('johndoe.csv', company='ftdna')

This command will write a file called 'genotype.txt' in the working directory; 
this contains your genotype in a format understood by DIYDodecad.

The 'standardize' command sorts SNP IDs alphabetically, removes 
headers/comments in the different company formats, and outputs a file 
'genotype.txt' where each line has four fields separated by whitespace:

SNPID CHROMOSOME POSITION GENOTYPE

6. 

You are now ready to launch DIYDodecad! You can do this from either your
operating system command prompt, or from within R.

a. At your operating system command prompt, go to the working directory (using 
the 'cd' command) and then enter:

In Windows:

DIYDodecadWin dv3.par

In 32bit Linux:

./DIYDodecadLinux32 dv3.par

In 64bit Linux:

./DIYDodecadLinux64 dv3.par

The program will now begin computing your admixture proportions for Dodecad v3 
components.

b. If you don't know how to use the operating system command prompt, you can 
launch DIYDodecad from within the R environment. Just type in at the R prompt:

system('DIYDodecadWin dv3.par')

'dv3' corresponds to Dodecad v3 results for Project participants. For more 
information:

https://spreadsheets.google.com/spreadsheet/ccc?key=0ArAJcY18g2GadDUyeEtjNnBmY09
EbnowN3M3UWRyNnc&hl=en_US&authkey=COCa89AJ
http://dodecad.blogspot.com/2011/06/design-of-dodecad-v3.html

The 2.0 version of the software allows for more detailed analysis of your 
genome by editing some options in the dv3.par file. See next section for 
details.

5) TECHNICAL DETAILS

DIYDodecad uses the EM algorithm to derive the maximum likelihood estimate of 
admixture proportions implementing the model of [1] as described by [2]. 

Each DIYDodecad calculator (e.g., 'dv3') consists of a series of files:

dv3.par		Parameters to run the program
dv3.alleles	Locus and allele information: SNPID MINOR_ALLELE MAJOR_ALLELE
dv3.12.F	Allele frequency data for putative ancestral populations
		These are in the format output by ADMIXTURE software
dv3.txt		Names of the ancestral populations

It is possible to make your own calculator, by creating and distributing a set 
of these four files. Note that the dv3.12.F and dv3.alleles files should be 
sorted alphabetically based on SNP IDs.

The parameter file (dv3.par for the 'dv3' calculator) includes the following 
lines:

1D-7
12
genotype.txt
166462
dv3.txt
dv3.12.F
dv3.alleles
verbose
genomewide

These are:

- Termination condition (maximum change of admixture proportions in successive 
iterations)
- Number of ancestral populations
- Input genotype file name
- Number of SNP markers
- List of ancestral population names
- Allele frequencies
- Allele information
- verbose/silent/progress
	silent: no intermediate output
	verbose: loglikelihood changes printed
	progress: intermediate solutions printed
- genomewide/bychr/byseg/target
	genomewide: admixture proportions across the entire genome
	bychr: admixture proportions separately for all chromosomes 1-22
	byseg: admixture proportions separately for regions within chromosomes
		NOTE: byseg mode requires two additional parameters
		window size:	how many SNPs a window consists of
		advance step: 	by how many SNPs to advance the window
	target: admixture proportions for a single region of the genome
		NOTE: target mode requires three additional parameters
		chromosome: 	1 to 22, the chromosome where the region is
		start:  	physical position (in bp) of the region's 
beginning
		end: 		physical position (in bp) of the region's end

Most of these you shouldn't alter, but you can change the following:

- Termination condition, e.g., to 1D-10 which will result in many more 
iterations run and potentially more accurate results, although in the vast 
majority of cases the results will be little affected
- verbose/silent/progress; in general, "verbose" is best. You might also try 
"progress" in "genomewide" mode, to see how the answer converges; in other 
modes iterations happen so fast (because of the smaller  number of SNPs) that 
make "progress" less useful
- genomewide/bychr/byseg/target. If you replace "genomewide" by "bychr", then 
your admixture proportions will be estimated per chromosome. If you replace 
"genomewide" with the following:

byseg
500
50

Then the program will use windows of 500 contiguous SNPs along a chromosome, 
and slide this window by increments of 50. So, it will first examine SNPs 1-500 
along chromosome 1, then 51-550, 101-600, and so on until the chromosome is 
exhausted, and will proceed to chromosome 2, and so on.

A small window size may result in more "noisy" results, as it becomes more 
difficult to estimate the origin of a segment of DNA of smaller size. On the 
other hand, it provides finer-scale resolution for segments of more remote 
ancestors. The program does not consider window sizes less than 10 SNPs.

A small advance step provides a more comprehensive scan of the genome, and 
results in many more segments being considered; the downside is that the 
running time of the analysis increases.

If you are interested in a particular genomic region, replace "genomewide" with

target
1
20000000
25000000

This will only examine the region on chromosome 1 between 20,000,000 and 
25,000,000 bp. Of course there may not be a SNP exactly at these positions, so 
the program will use the closest SNPs inside the window. In the above-given 
example, the target region will consist of 356 SNPs from 20008170 to 24953109bp.

The program will not consider regions with less than 10 SNPs; few SNPs make 
estimation of a segment's origin difficult, except for genomic regions with 
high inter-population frequency differences.

[1] H. Tang, J. Peng, P. Wang, N. Risch. Estimation of individual admixture: 
Analytical and study design considerations. Genet Epidemiol 28: 289-301, 2005.
[2] D.H. Alexander, J. Novembre, and K. Lange. Fast model-based estimation of 
ancestry in unrelated individuals. Genome Research, 19:1655-1664, 2009. 

6) RESULT FILES

Results are printed on screen and written in output files, depending on which  
mode was used (see Previous section):

genomewide.txt
bychr.txt
byseg.txt
target.txt

The bychr.txt and byseg.txt files are continually updated as the program runs, 
so you might want to open them and check out the results as they are being 
produced.

A note of caution should be kept in mind for the interpretation of small
admixture proportions: over a large region of your genome, such proportions may
hint at your having a smaller segment from a particular population hidden 
within that larger region; you may use a smaller window size to discover such
segments. 

Note, however, that as window sizes become smaller it is possible that some 
small admixture proportions can be noise, because there may not be enough SNPs
within the smaller region to differentiate sufficiently some of the 
populations.

7) VISUALIZATION

To aid in your exploration of your genome, I have included paint_byseg.r, a 
tool for visualizing admixture proportions along a chromosome. You can use this
in R by entering (in the working directory):

source('paint_byseg.r')

paint_byseg relies on the byseg.txt output file, so you should first carry byseg
analysis in DIYDodecad to create this file. 

You can use it by entering:

paint_byseg(chr=8,calc='dv3')

This will plot your Chromosome 8. It will probably be a very messy plot, so you
might want to limit yourself to a smaller region of the chromosome (in Mb), 
e.g.,

paint_byseg(chr=8, region=c(50,60), calc='dv3')

To make plots clearer, you may maximize their windows in R and use the Windows 
menu to switch between the command prompt and the plot. 

You can also write a plot to an image file:

paint_byseg(chr=8, region=c(50,60), calc='dv3', tofile=T)

or, to make it as big as you like:

paint_byseg(chr=8,region=c(50,60),calc='dv3', tofile=T, width=2000,height=1000)

This will output a file "Chromosome8.png" in your working directory.

To remove the clutter, you may choose to plot only the top/bottom few
components in a segment. You do this by:

paint_byseg(chr=8, region=c(50,60), calc='dv3', top=3)
paint_byseg(chr=8, region=c(50,60), calc='dv3', bottom=3)

Or, you may plot only the components you want, but be careful to spell them
correctly!

paint_byseg(chr=8, calc='dv3', choice=c("West_Asian", "East_European"))

Also, keep in mind that you can use the up and down arrows of your keyboard
to go back to your command history, so that you don't have to retype 
everything.

8) FEEDBACK

Direct all feedback about this program to dodecad@gmail.com. Bugs, suggestions 
for improvement, your results are all welcome.
