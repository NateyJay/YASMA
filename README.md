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
* This might be a good way to find miRNAs, but I think it would still be very problematic for fungi. Looking at the Tratr/Bocin annotations from Consuelo's work - virtually the entire genome is loci.


#### Expression thresholding
* Perhaps we could have a high minimum locus coverage threshold to find a locus, with less strict thresholds for it's edges. 
* I don't know if this would directly address the problems shortstack has with miRNA discovery.
* Binning (not this actually) or a density-driven algorithm could be a basis for this.
* A quantile-based RPM threshold?


#### miRNA-specific algorithms
* Highly-multiplexed folding would allow us to try to find resilient hairpins within loci. *Are these representative of real hairpins?* It looks like from the fungi hpRNA paper that we can recover miRNAs pretty well in control species using this approach. In any case, I think more folding is an answer to this problem.
* Maybe we could fold stretches that are highly stranded with multiple sizes. Most frequent hits could be condensed to find the most resilient hairpin.



## Current approach
  
The current version of the code is generally similar to ShortStack, but implements some of the changes above and some more suitable defaults to account for the changes. The pipeline currently follows 3 steps.
  
#### Assumptions
As an input for this pipeline, libraries should be trimmed to remove any trace of adapter sequences. As with all sRNA-seq analyses, untrimmed reads should be removed. Some size selection at the trimming step might be important, as reads shorter than ~15 nt become quite problematic for alignment. However, this tool relies on non-sRNA-related reads to also be included in the alignment to function properly.  
  
A good recommended trimming range would be retaining 15-30 nt reads, though this might need to be flexible depending on your organism of interest.

#### **1)** Alignment
Currently, this approach uses ShortStack as its aligner protocol (based on Bowtie1). This is most simply run with all (trimmed) libraries and the ```--alignment_only``` option. 

```
ShortStack --readfile wt_1.fa wt_2.fa mut_1.fa mut_2.fa \
--genomefile genome.fa \
--bowtie_cores 8 \
--alignment_only \
--cram \
--outdir alignment_dir
```


#### **2)** Precheck
This is meant to check the size distribution of reads, to help identify what are the likely sizes of dicer-derived small RNAs for the organism. In a well-known organim (like most plants or animals), you could likely just skip this step and provide a range for step 3 based on the literature.

```
annotate.py precheck -b ./alignment_dir/merged_alignments.cram \
-r wt_1 -r wt_2 \
-o annotation_complete
```

Produces:
```
  length	prop	prop_highest	abundance
  15	0.0	0.0001	169
  16	0.0	0.0001	128
  17	0.0	0.0001	107
  18	0.0001	0.0005	731
  19	0.018	0.0941	127,123
  20	0.0467	0.2452	331,033
  21	0.1218	0.639	862,850
  22	0.1907	1.0	1,350,239
  23	0.1169	0.613	827,728
  24	0.0735	0.3853	520,239
  25	0.0777	0.4074	550,146
  26	0.0738	0.3869	522,364
  27	0.0744	0.3904	527,094
  28	0.0702	0.368	496,933
  29	0.0701	0.3677	496,491
  30	0.0661	0.3466	467,936
```
Showing the sRNA length, proportion, and depth. This result clearly points to what we already knew from the literature: this species produces mostly 21 and 22 nt sRNAs. If this is more vague, you might use a more broad setting for the following step.
  

#### **3)** Annotation
Here the actual annotation takes place. Many options can affect the quality and discrimination of this step, so default and presets (a future option) are recommended for basic users. This approach basically takes 2 filters to define sRNA-producing regions and subsequently merges them to build loci.
  
Regions are defined by 2 rules:
1) Positions passing a minimal coverage of DCR-sized sRNAs (default 1.0 RPM). Regions with lower RPMs are also considered (default 1.0 RPM * 0.5), but only for extending loci which have contain a fully-passing region. 
2) Positions where a in a user-defined window (default 100 nt) there are more DCR-sized sRNAs overlapping by a certain fold (default 3 fold).
  
These regions are then merged based on a distance (default 150 nt), expanded to incorporate any overlapping reads, and padded to add some 'wiggle room' (default 10 nt).
  
  
To run, we must specify a few important options:  
```--bamfile (-b)``` for the merged alignment from before.  
  
```--output_directory (-o)``` to give a clear location for the output (default uses a UNIX-time ID).  
  
```--dicercall (-d)``` takes a list of sizes which you expect for dicer-derived sRNAs. This may be tricky with some samples. Inputs here will interpolate a list (you only need to provide the start and finish).  
  
```--readgroups (-r)``` takes a list of the read-groups (libraries) from which you want to form the annotation. For example, you might want to exclude libraries which have a sRNA-biogenesis mutation, otherwise you would would specifically annotate non-sRNAs. Read-groups excluded here will still be included in the counts table. 
  
  
Here's an example call:
```
annotate.py annotate -b ./alignment_dir/merged_alignments.cram \
-d 21 -d 22 \
-r wt_1 -r wt_2 \
-o annotation_complete
```
Will produce a full annotation in the output folder.  
  
  
A full description of the options is as follows:  
```
annotate.py annotate --help
Usage: annotate.py annotate [OPTIONS]

  Main annotation suite.

Options:
  -a, --alignment_file PATH       Alignment file input (bam or cram).
                                  [required]
  -r, --annotation_readgroups TEXT
                                  List of read groups (RGs, libraries) to be
                                  considered for the annotation. 'ALL' uses
                                  all readgroups for annotation, but often
                                  pertainent RGs will need to be specified
                                  individually.  [required]
  -d, --dicercall INTEGER         List of sRNA lengths that derive from dicer.
  -o, --output_directory PATH     Directory name for annotation output
  -f, --force                     Force remake of supporting files
  --partial_wigs                  Only make wiggle files associated with
                                  essential functions (ignoring size and
                                  strand specific coverages. (May improve
                                  speed)
  --window INTEGER                Window size (centered on position) for
                                  determining DCR vs non-DCR read ratio
                                  (counting overlapping reads).
  --merge_dist INTEGER            Maximum gap size between valid regions to
                                  merge to a single locus.
  --pad INTEGER                   Number of bases arbitrarily added to either
                                  end of a defined locus.
  --rpm_cutoff FLOAT              RPM depth threshold for DCR-sized reads to
                                  be considered as a valid region.
  --extension_ratio FLOAT         Fraction of RPM threshold to be considered
                                  for extending a locus boundaries
  --dicer_ratio FLOAT             Ratio of dicer to non-dicer reads to be
                                  considered for a valid region
  --help                          Show this message and exit.


```

### Output files
The annotator produces several key output files:


```Log.txt``` is a text file logging all of the runtime messages.  
  
```Results.txt``` is a tab-delimited description of all annotated loci.  
  
```Counts.txt``` is a tab-delimited table of the read depth for each read group in each locus. Counts only consider reads entirely enclosed in a locus.  
  
```Annotation.gff``` is a gff-formatted annotation of all loci. If prerequisites are installed, the tool will automatically sort, zip, and index the annotation for jbrowse use.  
  
```TopReads.txt``` is a tab-delimited table of the most abundant reads for all loci. These are given as ranked sequences, showing their depth, rpm, and cumulative proportion of a locus depth. 
  
```wigs``` and ```bigwigs``` are folders containing wig and bigwig (assuming prerequisites) coverage files. This includes coverages for all DCR sizes by strand, DCR and non-DCR coverages (ignoring strand), and metrics for passing and considering a sRNA-regions.

  
## Planned features

* **Locus input.** In addition to discovery, this tool will take input loci which can be evaluated independently (or along side) *de novo* annotation.
* **miRNA discovery methods.** This will focus on folding stranded regions of loci, including outside of locus boundaries, to try to find stable hairpins which meet miRNA rules. These will be reported as separate, possibly "daughter" loci.
* **Preset settings.** These can be used to mimic different annotation rules - for example shortstack or the default settings for this software.
* **Dicer agnositic settings.** This requires more thinking, as dicer-dependency is an essential aspect of this tool. However, there will likely be examples where it is hard or impossible to discern a proper dicer-range. It may be important to have informed settings for this tool to handle annotation lacking dicer information.
* **JBrowse configuration.** Automatically producing configuration file for jbrowse.





