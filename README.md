# annotator
 

A general look at better ways to define loci in sRNA seq data.


### Problems with current approaches
* ShortStack does a good job, with a easily understood algorithm and set of rules. However, it misses lots of miRNAs because of concatenation of loci.
* This problem is exacerbated in fungi - high background sRNA levels make it hard to distinguish off-sized, degradation-related, and other loci from dicer or RNAi-related loci.
* Moreover, ShortStack defines the boundaries of loci based on a merge window that require only a single read to extend. This means in libraries which are overly-deep compared to the genome size, loci can be come trivially large (annotating the whole genome). This is what we see in several fungi.
* **The problems we see in fungi should be the emphasis of this software**


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

The current version of the code is generally similar to ShortStack, but implements some of the changes above and some more suitable defaults to account for the changes. The pipeline currently follows 3 steps: 

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
01-precheck.py -b ./alignment_dir/merged_alignments.cram -o precheck
```

Produces:
```
  15	0.005	622,756
  16	0.005	595,494
  17	0.007	808,679
  18	0.009	1,053,851
  19	0.022	2,665,887
  20	0.043	5,238,861
  21	0.115	14,118,431
  22	0.158	19,344,300
  23	0.079	9,648,308
  24	0.071	8,717,748
  25	0.07	8,622,328
  26	0.063	7,677,058
  27	0.064	7,867,814
  28	0.063	7,769,189
  29	0.064	7,800,495
  30	0.056	6,932,318
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
02-annotate.py -b ./alignment_dir/merged_alignments.cram \
-d 21 22 \
-r wt_1 wt_2 \
-o annotation_complete
```
Will produce a full annotation in the output folder.  
  
A full description of the options is as follows:  
```
usage: 02-annotate.py [-h] -b [BAM_FILE] [-o [OUTPUT_DIRECTORY]] [-d DICERCALL [DICERCALL ...]] [-r READGROUPS [READGROUPS ...]] [-f] [--partial_wigs] [--window [WINDOW]] [--merge_dist [MERGE_DIST]]
                      [--pad [PAD]] [--rpm [RPM]] [--dicer_ratio [DICER_RATIO]] [--extension_ratio [EXTENSION_RATIO]]

options:
  -h, --help            show this help message and exit
  -b [BAM_FILE], --bam_file [BAM_FILE]
                        bamfile of aligned small RNAs (tested with shortstack)
  -o [OUTPUT_DIRECTORY], --output_directory [OUTPUT_DIRECTORY]
                        working folder for locus analysis
  -d DICERCALL [DICERCALL ...], --dicercall DICERCALL [DICERCALL ...]
                        list of sRNA sizes to be analyzed separately
  -r READGROUPS [READGROUPS ...], --readgroups READGROUPS [READGROUPS ...]
                        list of read groups (libraries) to be considered for the annotation. All present read groups will be considered for counting.
  -f, --force           force remake of supporting files
  --partial_wigs        only make wiggle files associated with essential functions (ignoring size and strand specific coverages)
  --window [WINDOW]     kernel desity bandwith window
  --merge_dist [MERGE_DIST]
                        maximum gap size between valid regions to merge to a single locus
  --pad [PAD]           number of bases added to either end of a defined locus (arbitrary)
  --rpm [RPM]           rpm depth threshold for nucleating a locus
  --dicer_ratio [DICER_RATIO]
                        ratio of dicer to non-dicer reads to be considered for a locus region
  --extension_ratio [EXTENSION_RATIO]
                        fraction of rpm threshold to be considered for extending a locus' boundaries

```

### Output files
The annotator produces several key output files:


```Log.txt``` is a text file logging all of the runtime messages.  
  
```Results.txt``` is a tab-delimited description of all annotated loci.  
  
```Counts.txt``` is a tab-delimited table of the read depth for each read group in each locus. Counts only consider reads entirely enclosed in a locus.  
  
```Annotation.gff``` is a gff-formatted annotation of all loci. If prerequisites are installed, the tool will automatically sort, zip, and index the annotation for jbrowse use.  
  
```TopReads.txt``` is a tab-delimited table of the most abundant reads for all loci. These are given as ranked sequences, showing their depth, rpm, and cumulative proportion of a locus depth. 
  
```Coverages``` is a folder containing wig and bigwig (assuming prerequisites) coverage files. This includes coverages for all DCR sizes and strands; DCR and non-DCR coverages grouped and non-stranded; metrics for passing and considering a sRNA-region.









