# Day 3 Morning
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On day 1, we ran through a pipeline to map reads against a reference genome and call variants, but didn’t do much with the variants we identified. Among the most common analyses to perform on a set of variants is to construct phylogenetic trees. Here we will explore different tools for generating and visualizing phylogenetic trees, and also see how recombination can distort phylogenetic signal.

For the first several exercises, we will use the A. baumannii genomes that we worked with yesterday afternoon. 
The backstory on these genomes is that Abau_A, Abau_B and Abau_C are representatives of three clones (as defined by pulsed-field gel electrophoresis - a low-resolution typing method) that were circulating in our hospital. 

One of the goals of our published study was to understand the relationship among these clones to discern whether: 

1) the three clones represent three independent introductions into the hospital or 

2) the three clones originated from a single introduction into the hospital, with subsequent genomic rearrangement leading to the appearance of unique clones. 

The types of phylogenetic analyses you will be performing here are the same types that we used to decipher this mystery.
The other two genomes you will be using are ACICU and AB0057. ACICU is an isolate from a hospital in France, and its close relationship to our isolates makes it a good reference for comparison. AB0057 is a more distantly related isolate that we will utilize as an out-group in our phylogenetic analysis. The utility of an out-group is to help us root our phylogenetic tree, and gain a more nuanced understanding of the relationship among strains.

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```
wd

# or

cd /scratch/micro612w17_fluxod/username

cp -r /scratch/micro612w17_fluxod/shared/data/day3_morn ./

```

## Perform whole genome alignment with [Mauve](http://darlinglab.org/mauve/mauve.html) and convert alignment to other useful formats
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

An alternative approach for identification of variants among genomes is to perform whole genome alignments of assemblies. If the original short read data is unavailable, this might be the only approach available to you. Typically, these programs don’t scale well to large numbers of genomes (e.g. > 100), but they are worth being familiar with. We will use the tool mauve for constructing whole genome alignments of our five A. baumannii genomes.

>i. Perform mauve alignment and transfer xmfa back to flux

Use sftp to get genomes onto your laptop

```
Run these commands on your local system/terminal:

cd ~/Desktop (or wherever your desktop is) 

mkdir Abau_mauve

cd Abau_mauve 

> Now copy Abau_genomes folder residing in your day3_morn folder using scp or sftp:

scp -r username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/Abau_genomes ./

OR

sftp –r username@flux-login.arc-ts.umich.edu 
cd /scratch/micro612w17_fluxod/username/day3_morn 
get Abau_genomes

```

Run mauve to create multiple alignment

```

i. Open mauve 
ii. File -> align with progressiveMauve 
iii. Click on “Add Sequnce” and add each of the 5 genomes you just downloaded
iv. Name the output file “mauve_ECII_outgroup” and make sure it is in the directory you created for this exercise 
v. Click Align! 
vi. Wait for Mauve to finish and explore the graphical interface

```

Use sftp or scp to transfer your alignment back to flux for some processing

```

cd ~/Desktop/Abau_mauve
sftp –r username@flux-login.arc-ts.umich.edu 
cd /scratch/micro612w17_fluxod/username/day3_morn 
put mauve_ECII_outgroup

OR

scp ~/Desktop/Abau_mauve/mauve_ECII_outgroup username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn 

```
 
>ii. Convert alignment to fasta format

Mauve produces alignments in .xmfa format (use less to see what this looks like), which is not compatible with other programs we want to use. We will use a custom script convert_msa_format.pl to change the alignment format to fasta format

<!-- correction pending-->
```

Now run these command in day3_morn folder on flux:

module load bioperl

perl convert_msa_format.pl -i mauve_ECII_outgroup -o mauve_ECII_outgroup.fasta -f fasta -c

```

## Perform some DNA sequence comparisons and phylogenetic analysis in [APE](http://ape-package.ird.fr/), an R package
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

There are lots of options for phylogenetic analysis. Here, we will use the ape package in R to look at our multiple alignments and construct a tree using the Neighbor Joining method. 

Note that ape has a ton of useful functions for more sophisticated phylogenetic analyses!

>i. Get fasta alignment you just converted to your own computer using sftp or scp

```

cd ~/Desktop/Abau_mauve

sftp –r username@flux-login.arc-ts.umich.edu 
cd /scratch/micro612w17_fluxod/username/day3_morn 
get mauve_ECII_outgroup.fasta

OR

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/mauve_ECII_outgroup.fasta ./

```

ii. Read alignment into R

Fire up RStudio and install/load ape

Use the read.dna function in ape to read in you multiple alignments. 
Print out the variable to get a summary.

```
install.packages("ape")
library(ape)
abau_msa = read.dna('mauve_ECII_outgroup.fasta', format = "fasta") 
```

>iii. Get variable positions

The DNA object created by read.dna can also be addressed as a matrix, where the columns are positions in the alignment and rows are your sequences. We will next treat our alignment as a matrix, and use apply and colSums to get positions in the alignment that vary among our sequences. Examine these commands in detail to understand how they are working together to give you a logical vector indicating which positions vary in your alignment.

```

abau_msa_bin = apply(abau_msa, 2, FUN = function(x){x == x[1]}) 

abau_var_pos = colSums(abau_msa_bin) < 5
```

>iv. Get non-gap positions

For our phylogenetic analysis we want to focus on the core genome, so we will next identify positions in the alignment where all our genomes have sequence.

```
non_gap_pos = colSums(as.character(abau_msa) == '-') == 0
```

>v. Count number of variants between sequences

Now that we know which positions in the alignment are core and variable, we can extract these positions and count how many variants there are among our genomes. Do count pairwise variants we will use the dist.dna function in ape. The model parameter indicates that we want to compare sequences by counting differences. Print out the resulting matrix to see how different our genomes are.

```

abau_msa_var = abau_msa[,abau_var_pos & non_gap_pos ]
var_count_matrix = dist.dna(abau_msa_var, model = "N")

```

>vi. Construct phylogenetic tree

Now we are ready to construct our first phylogenetic tree! 

We are going to use the Neighbor Joining algorithm, which takes a matrix of pairwise distances among the input sequences and produces the tree with the minimal total distance. In essence, you can think of this as a distance-based maximum parsimony algorithm, with the advantage being that it runs way faster than if you were to apply a standard maximum parsimony phylogenetic reconstruction.

As a first step we are going to build a more accurate distance matrix, where instead of counting variants, we will measure nucleotide distance using the Jukes-Cantor model of sequence evolution. This is the simplest model of sequence evolution, with a single mutation rate assumed for all types of nucleotide changes.

```
dna_dist_JC = dist.dna(abau_msa, model = "JC")
```

Next, we will use the ape function nj to build our tree from the distance matrix

```
abau_nj_tree = nj(dna_dist_JC)
```

Finally, plot your tree to see how the genomes group.

```
plot(abau_nj_tree)
```

## Perform SNP density analysis to discern evidence of recombination
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

An often-overlooked aspect of a proper phylogenetic analysis is to exclude recombinant sequences. Homologous recombination in bacterial genomes is a mode of horizontal transfer, wherein genomic DNA is taken up and swapped in for a homologous sequence. The reason it is critical to account for these recombinant regions is that these horizontally acquired sequences do not represent the phylogenetic history of the strain of interest, but rather in contains information regarding the strain in which the sequence was acquired from. One simple approach for detecting the presence of recombination is to look at the density of variants across a genome. The existence of unusually high or low densities of variants is suggestive that these regions of aberrant density were horizontally acquired. Here we will look at our closely related A. baumannii genomes to see if there is evidence of aberrant variant densities.

>i. Subset sequences to exclude the out-group

For this analysis we want to exclude the out-group, because we are interested in determining whether recombination would hamper our ability to reconstruct the phylogenetic relationship among our closely related set of genomes.  

>Note that the names of the sequences might be different for you, so check that if the command doesn’t work.

```

abau_msa_no_outgroup = abau_msa[c('ACICU_genome','AbauA_genome','AbauC_genome','AbauB_genome'),]

```

>ii. Get variable positions

Next, we will get the variable positions, as before

```

abau_msa_no_outgroup_bin = apply(abau_msa_no_outgroup, 2, FUN = function(x){x == x[1]}) 

abau_no_outgroup_var_pos = colSums(abau_msa_no_outgroup_bin) < 4

```

>iii. Get non-gap positions

Next, we will get the core positions, as before

```

abau_no_outgroup_non_gap_pos = colSums(as.character(abau_msa_no_outgroup) == '-') == 0

```

>iv. Create overall histogram of SNP density

Finally, create a histogram of SNP density across the genome. Does the density look even, or do you think there might be just a touch of recombination?

```
hist(which(abau_no_outgroup_var_pos & abau_no_outgroup_non_gap_pos), 10000)
```

## Perform recombination filtering with [Gubbins](https://www.google.com/search?q=gubbins+sanger&ie=utf-8&oe=utf-8)
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

Now that we know there is recombination, we know that we need to filter out the recombinant regions to discern the true phylogenetic relationship among our strains. In fact, this is such an extreme case (~99% of variants of recombinant), that we could be totally misled without filtering recombinant regions. To accomplish this we will use the tool gubbins, which essentially relies on elevated regions of variant density to perform recombination filtering.

>i. Run gubbins on your fasta alignment

Go back on flux and load modules required by gubbins
<!-- Correction pending -->

```
Check if gubbins run after loading newer version flux modules

Older version:
module load python/2.7.3 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins

Newer version:
module load python-anaconda2/201607 biopython dendropy reportlab fasttree RAxML fastml/gub gubbins

```

Run gubbins on your fasta formatted alignment

```
d3m

# or

cd /scratch/micro612w17_fluxod/username/day3_morn

run_gubbins.py -v -f 50 -o Abau_AB0057_genome mauve_ECII_outgroup.fasta

```

>ii. Create gubbins output figure

Gubbins produces a series of output files, some of which can be run through another program to produce a visual display of filtered recombinant regions. Run the gubbins_drawer.py script to create a pdf visualization of recombinant regions. 

The inputs are: 

1) the recombination filtered tree created by gubbins (mauve_ECII_outgroup.final_tree.tre),

2) the pdf file to create (mauve_ECII_outgroup.recombination.pdf) and 

3) a .embl representation of recombinant regions (mauve_ECII_outgroup.recombination_predictions.embl).

```

gubbins_drawer.py -t mauve_ECII_outgroup.final_tree.tre -o mauve_ECII_outgroup.recombination.pdf mauve_ECII_outgroup.recombination_predictions.embl

```
>iii. Download and view gubbins figure and filtered tree

Use sftp or scp to get gubbins output files into Abau_mauve on your local system

```

cd ~/Desktop/Abau_mauve

sftp –r username@flux-login.arc-ts.umich.edu 
cd /scratch/micro612w17_fluxod/username/day3_morn 
get mauve_ECII_outgroup.recombination.pdf 
get mauve_ECII_outgroup.final_tree.tre

OR

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/mauve_ECII_outgroup.recombination.pdf  ./
scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/mauve_ECII_outgroup.final_tree.tre  ./

```

Open up the pdf and observe the recombinant regions filtered out by gubbins. Does it roughly match your expectations based upon your SNP density plots?

Finally, lets look at the recombination-filtered tree to see if this alters our conclusions. 

To view the tree we will use [Seaview](http://doua.prabi.fr/software/seaview), which is a multi-purpose tool for: 

1) visualization/construction of multiple alignments and 

2) phylogenetic tree construction. 

Here, we will just use Seaview to view our gubbins tree.

```

In seaview: 

Go to Trees -> import tree (mauve_ECII_outgroup.final_tree.tre) 
To view sub-tree of interest click on “sub-tree” and select the sub-tree excluding the out-group

```


How does the structure look different than the unfiltered tree?

> Note that turning back to the backstory of these isolates, Abau_B and Abau_C were both isolated first from the same patient. So this analysis supports that patient having imported both strains, which likely diverged at a prior hospital at which they resided.

## Create annotated publication quality trees with [iTOL](http://itol.embl.de/)
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

For the final exercise we will use a different dataset, composed of USA300 methicillin-resistant Staphylococcus aureus genomes. USA300 is a strain of growing concern, as it has been observed to cause infections in both hospitals and in otherwise healthy individuals in the community. An open question is whether there are sub-clades of USA300 in the hospital and the community, or if they are all the same. Here you will create an annotated phylogenetic tree of strains from the community and the hospital, to discern if these form distinct clusters.

>i. Download MRSA genome alignment from flux

Use sftp or scp to get genomes onto your laptop

```

cd ~/Desktop (or wherever your desktop is) 
mkdir MRSA_genomes 
cd MRSA_genomes

sftp –r username@flux-login.arc-ts.umich.edu 
cd /scratch/micro612w17_fluxod/username/day3_morn 
get 2016-3-9_KP_BSI_USA300.fa 
get 2016-3-9_KP_BSI_USA300_iTOL_HA_vs_CA.txt

OR

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/2016-3-9_KP_BSI_USA300.fa  ./
scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day3_morn/2016-3-9_KP_BSI_USA300_iTOL_HA_vs_CA.txt  ./


```

>ii. Look at SNP density for MRSA alignment in R

Before we embark on our phylogenetic analysis, lets look at the SNP density to verify that there is no recombination

```

mrsa_msa = read.dna('2016-3-9_KP_BSI_USA300.fa', format = 'fasta') 
mrsa_msa_bin = apply(mrsa_msa, 2, FUN = function(x){x == x[1]}) 
mrsa_var_pos = colSums(mrsa_msa_bin) < nrow(mrsa_msa_bin) 
hist(which(mrsa_var_pos), 10000)

```

Does it look like there is evidence of recombination?

>iii. Create fasta alignment with only variable positions

Next, lets create a new fasta alignment file containing only the variant positions, as this will be easier to deal with in Seaview

```

write.dna(mrsa_msa[, mrsa_var_pos], file = '2016-3-9_KP_BSI_USA300_var_pos.fa', format = 'fasta')

```

>iv. Read alignment into Seaview and construct Neighbor Joining tree

In the previous exercise, we used Seaview to look at a pre-existing tree, here we will use Seaview to create a tree from a
multiple sequence alignment 

Read in multiple alignment of variable positions

```
Go to File -> open ('2016-3-9_KP_BSI_USA300_var_pos.fa)
```

Construct Neighbor Joining phylogenetic tree with default parameters (note, this will take a few minutes)

```
Go to Trees -> select Distance Methods -> BioNJ -> (Select Bootstrap with 20 replicates) -> Go
```

Save your tree

```
File -> Save rooted tree
```

Note that in your research it is not a good idea to use these phylogenetic tools completely blind and I strongly encourage embarking on deeper learning yourself, or consulting with an expert before doing an analysis for a publication

v. Read tree into iTOL

```

To make a prettier tree and add annotations we will use iTOL (http://itol.embl.de/). 

Go to http://itol.embl.de/

To load your tree, click on upload, and select the rooted tree you just created in Seaview

```

Explore different visualization options for your tree (e.g. make it circular, show bootstrap values, try collapsing nodes/branches)

Note that you can always reset your tree if you are unhappy with the changes you’ve made

>vi. Add annotations to tree

One of the most powerful features of iTOL is its ability to overlay diverse types of descriptive meta-data on your tree (http://itol.embl.de/help.cgi#datasets). Here, we will overlay our data on whether an isolate was from a community or hospital infection. To do this simply drag-and-drop the annotation file (2016-3-9_KP_BSI_USA300_iTOL_HA_vs_CA.txt) on your tree and voila! 

> Do community and hospital isolates cluster together, or are they inter-mixed?

