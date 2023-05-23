# YASMA
*<ins>Y</ins>et <ins>a</ins>nother <ins>sm</ins>all RNA <ins>a</ins>nnotator*
 
**warning**: this project is in active development and the readme is woefully out of date. This tool follows a much different syntax in the current version.

### Problems with current approaches
* ShortStack does a good job, with a easily understood algorithm and set of rules. However, it misses lots of miRNAs because of concatenation of loci.  
  
* This problem is exacerbated in fungi - high background sRNA levels make it hard to distinguish off-sized, degradation-related, and other loci from dicer or RNAi-related loci.  
  
* Moreover, ShortStack defines the boundaries of loci based on a merge window that require only a single read to extend. This means in libraries which are overly-deep compared to the genome size, loci can be come trivially large (annotating the whole genome). **This is what we see in several fungi.**



### Ideas for implementation
*Some implemented, others not...*

#### Loci and sub-loci.
* one approach could be to try to only provide more detail to those defined by shortstack. Here, we would define loci which are within shortstack loci, hopefully finding shorter loci which are miRNAs.
* This might be a good way to find miRNAs, but I think it would still be very problematic for fungi. Looking at the Tratr/Bocin annotations - virtually the entire genome is loci.


#### Expression thresholding
* Perhaps we could have a high minimum locus coverage threshold to find a locus, with less strict thresholds for it's edges.
* I don't know if this would directly address the problems shortstack has with miRNA discovery.
* Binning (not this actually) or a density-driven algorithm could be a basis for this.
* A quantile-based RPM threshold?
* A poisson-probability-based threshold.


#### miRNA-specific algorithms
* Highly-multiplexed folding would allow us to try to find resilient hairpins within loci. *Are these representative of real hairpins?* It looks like from the fungi hpRNA paper that we can recover miRNAs pretty well in control species using this approach. In any case, I think more folding is an answer to this problem.
* Maybe we could fold stretches that are highly stranded with multiple sizes. Most frequent hits could be condensed to find the most resilient hairpin.


# YASMA functionality

The tool is divided into 3 general categories of functions. 

* **Annotation** provides the main backbone of annotation for the tool. Currently, there are 2 approaches: 
1) using a poisson distribution to determine thresholding, and 
2) using a priori description of dicer-derived-sRNA sizes to pick regions which produce significantly more sRNAs.

The poisson approach appears to work the best (by far) and will be used for further testing.

* **Calculation** contains several tools to produce a fully usable annotation.
* **Utilities** provide some convenience functions for the user.



## Main (poisson) annotation pipeline
  
The current version of the code is generally similar to ShortStack, but implements some of the changes above and some more suitable defaults to account for the changes. The pipeline currently follows 3 steps.
  
#### Assumptions

#### **1)** Pre-processing

As an input for this pipeline, libraries should be trimmed to remove any trace of adapter sequences. As with all sRNA-seq analyses, untrimmed reads should be removed. Some size selection at the trimming step might be important, as reads shorter than ~15 nt become quite problematic for alignment. This tool also relies on non-sRNA-related reads (usually short or long) to also be included in the alignment to function properly.
  
A good recommended trimming range would be retaining 15-30 nt reads, though this might need to be flexible depending on your organism of interest.

I usually use the tool cutadapt with the following options:
```
cutadapt \
-a [adapter_seq] \
--minimum-length 10 \
-O 4 \
--discard-untrimmed \
--max-n 0 \
-o [trimmed_out_file] \
[untrimmed_file]
```


#### **2)** Alignment
Currently, this approach uses ShortStack as its aligner protocol (based on Bowtie1). This is most simply run with all (trimmed) libraries and the ```--alignment_only``` option. *Note: this command is written for ShortStack 3.X*

```
ShortStack --readfile wt_1.fa wt_2.fa mut_1.fa mut_2.fa \
--genomefile genome.fa \
--bowtie_cores 8 \
--alignment_only \
--cram \
--outdir alignment_dir
```



#### **3)** Annotation
Here the actual annotation takes place. Many options can affect the quality and discrimination of this step, so default and presets (a future option) are recommended for basic users. 



**Calculating poisson** model by chromosome

Poisson-distributions are well suited to sRNA alignment analyses, as alignment counts are integers and occur with frequency over a genomic space.

This model first scans through to alignment for each chromosome, calculating the lambda coefficient. To reflect this most accurately, the tool uses gene-annotations as an input, as these will produce lots of sRNA-sized reads as the result of degradation. These can produce high background in fungi, so YASMA focuses on intergenic regions (as can be most easily defined).

**Scoring window probabilities**

Next, YASMA iterates through the alignment, and scores windows (default 40 nt) based on their poisson probability in a given chromosome. By default, it considers regions with P < 0.00001 as candidates. Candidate regions are then expanded until they find a window where depth is smaller than lambda - this is the preliminary boundary.

**Merging and trimming**

Passing regions are then merged if they are less than a threshold distance apart (default 150 nt). Finally, merged regions are trimmed to have depths at their edges with are not lower than 0.05 * their peak depth. This helps make highly peaky loci more concise.


**Repeat for all chromosomes**

This process is then repeated for every scaffold/contig/chromosome in the alignment.

### Output files and post processing









<!-- 
### Output files
The annotator produces several key output files:


```Log.txt``` is a text file logging all of the runtime messages.  
  
```Results.txt``` is a tab-delimited description of all annotated loci.  
  
```Counts.txt``` is a tab-delimited table of the read depth for each read group in each locus. Counts only consider reads entirely enclosed in a locus.  
  
```Annotation.gff``` is a gff-formatted annotation of all loci. If prerequisites are installed, the tool will automatically sort, zip, and index the annotation for jbrowse use.  
  
```TopReads.txt``` is a tab-delimited table of the most abundant reads for all loci. These are given as ranked sequences, showing their depth, rpm, and cumulative proportion of a locus depth. 
  
```wigs``` and ```bigwigs``` are folders containing wig and bigwig (assuming prerequisites) coverage files. This includes coverages for all DCR sizes by strand, DCR and non-DCR coverages (ignoring strand), and metrics for passing and considering a sRNA-regions. -->

  
## Planned features

* **Locus input.** In addition to discovery, this tool will take input loci which can be evaluated independently (or along side) *de novo* annotation.
* **miRNA discovery methods.** This will focus on folding stranded regions of loci, including outside of locus boundaries, to try to find stable hairpins which meet miRNA rules. These will be reported as separate, possibly "daughter" loci.
* **Preset settings.** These can be used to mimic different annotation rules - for example shortstack or the default settings for this software.
* **Dicer agnositic settings.** This requires more thinking, as dicer-dependency is an essential aspect of this tool. However, there will likely be examples where it is hard or impossible to discern a proper dicer-range. It may be important to have informed settings for this tool to handle annotation lacking dicer information.
* **JBrowse configuration.** Automatically producing configuration file for jbrowse.





