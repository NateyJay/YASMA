
<p align="center"><img src="images/logo2.png" alt="main logo" width="300" /></p>

YASMA is a sRNA annotation suite based around the tool **YASMA-tradeoff**. This tool performs **genome-based** *de novo* annotations which is focused on the ***total small RNA population***, not just a few classes like miRNAs.

In addition to the annotator, YASMA has many helper functions to facilitate reproducible analysis of sRNAs. These support the annotator, but are not required, as YASMA-tradeoff can produce annotations from any small RNA alignment (.bam).

### Published work available

This work is now [published](https://doi.org/10.1016/j.csbj.2025.05.045)! Please have a look at the publication if you want to know more about the benefits of YASMA-tradeoff and the lessons we learned about sRNAs in fungi.

## Overview

[Background](#background)  
[Installation](doc/installation.md)  
[Quick start](#quickstart)  
[YASMA basics](#basics)  
[Modules](#modules)  
[Cookbook](#cookbook)  


<a name="background"/>

## Background

#### Other tools for this?

There are other approaches that follow a similar strategy, namely [ShortStack](https://github.com/MikeAxtell/ShortStack). Yasma was born out of solving some persistent issues with ShortStack3, which appear to be exacerbated in challenging systems, such as sRNAs in Fungi (ShortStack4 has resolved some of these).



#### Problems with current approaches
* **Creeping annotations** which don't don't represent the shape of a expressed region, sometimes leading to very large loci.
* **Under-merging** where numerous similar loci are annotated separately due to sequencing gaps.
* **Weak identifications** of sRNA locus classes in non-canonically sized loci and organisms where we know relatively little (fungi, for example).



#### Balancing sensitivity and specificity


Annotation based on these alignments follows a multi step approach:
1. Building an average ***sRNA coverage profile*** from replicate libraries.
2. Calculating the RPM threshold which best balances annotating the most reads in the smallest genomic space.
3. Building a profile of genomic-regions which are sufficiently deep based on this threshold.
4. Merging of peaks which have ***similar sRNA profiles***.
5. Filtering loci to remove those which might be real, but don't have sufficient depth, density, or complexity to assess them.

This results in contiguous loci which are more homogenous in profile. It also tends to avoid over-annotation of background sequences.

</a>

<a name="quickstart"/>


## Quick start

YASMA is now available on bioconda! This makes installing the tool and dependencies quite streamlined.

```
conda create -n my_env
conda activate my_env
conda install bioconda::yasma

yasma -h
```


Yasma is organized into several modules, made with the CLI-module [click](https://click.palletsprojects.com/). These modules are organized into several major sections which are generally ordered by processing step:

```
Options:
  -h, --help  Show this message and exit.

Commands:

  Preliminary:
    inputs        Initialize a project and log inputs for later analyses

  Processing:
    download      Download libraries from the NCBI SRA using their SRR code
    adapter       Tool to check untrimmed-libraries for 3' adapter content.
    trim          Wrapper for trimming using cutadapt.
    align         Aligner based on shortstack3-style weighting

  Annotation:
    tradeoff      Annotator using focused capturing the most reads in the...

  Calculation:
    context       Compares annotations to identify cluster genomic context.
    count         Gets counts for all readgroups, loci, strand, and sizes.
    hairpin       Evaluates annotated loci for hairpin or miRNA structures.
    jbrowse       Build coverage and config files for jbrowse2
    coverage      Produces combined bigwig coverage files

  Utilities:
    size-profile  Convenience function for calculating aligned size profile.
```

An absolutely simple description of how to run this analysis is like follows, using hypothetical libraries in the SRA. *This assumes you have added `yasma.py` to your PATH variable*

```
yasma.py inputs -o my_analysis -s SRR1111111 SRR2222222 -g /path/to/my/genome.fa

cd ./my_analysis

# preprocessing
yasma.py download
yasma.py adapter
yasma.py trim

# alignment
yasma.py align

# annotation
yasma.py tradeoff

# reporting on annotation
yasma.py count
yasma.py hairpin
yasma.py jbrowse -j /path/to/jbrowse2/directory/
```

</a>

<a name="basics"/>

## YASMA basics
### Directory oriented analysis

To help with ease of use, Yasma orients all of its analyses around a directory. Files produced and referenced by yasma are all stored in the `inputs.json` file, using relative paths. Analyses that produce outputs will automatically update this file, meaning you need not manually transmit information from one module to the next (for example: finding an adapter sequence, then trimming the libraries with it). 

`inputs.json` is human-readable and can be pretty easily modified manually, though not normally advisable.

This means that you need to initialize a directory, by specifiying an `--output_directory`/`-o`. Once this is done, you can run yasma from this directory without a problem. Any Yasma module can initialize (making an `inputs.json`), so you can skip to tradeoff if you have an alignment ready. You can even specify the current directory with `.`, which will import this directory's name as the project name.

Modules can be run simply with `yasma.py [module] [...]`, and modules should tell you what inputs or requirements you are lacking.


### Preliminary step - *inputs*

Inputs is not a required step, but it can be a major time-saver. Basically, it produces a file `config.json` which can store all input files for your analysis. If you specify these files here, you need not call them in subsequent steps. 

This also lets you know if there are incongruities in your data. For example, it compares chromosome names found in your reference genome to a gene annotation, showing if they don't match. This can frequently be a real time-saver as it catches common errors.



</a>

<a name="modules"/>

## Modules
### Processing

Using these modules, yasma can look for adapter sequences, trim libraries (using cutadapt), and align them to a genome. Alignment mimics the process of ShortStack3/4, utilizing bowtie1 as the alignment engine and placing multimapping reads based on unique-mapper weighting. Other alignments (bam format) may be used, but the `@RG` field is required in the bam (natively produced in YASMA or ShortStack alignments).


### Annotation

The main annotation module is called `tradeoff`, due to its threshold finding with a read vs genome tradeoff. This analysis should work on any shortstack bam/cram alignment. 

There are several options, but an essential one specified here is `-r, --annotation_readgroups`. This allows you to make your annotation based on a smaller group of libraries from your whole alignment. This is really useful when working with a large analysis with multiple replicate-groups, and you might only want to annotate sRNAS in one of them (e.g. wt replicates among many mutants).

Many outputs are produced from this step, with some described here:

* `loci.gff3` and `loci.txt` - the core annotation output, identifying loci and their dimensions (in gff and tabular formats). 
* `removed_loci.gff3` - loci that have been filtered in one of the locus filter steps.
* `coverage.bw` and `kernel.bw` - track files associated with the sRNA alignment coverage and padded_coverage by max() (kernel).
* `regions.gff3` and `revised_regions.gff3` - annotation files of distinct sRNA regions based off the padded_coverage, and the revision of those regions including nearby similar sRNAs.
* `thresholds.txt` - a table showing the percent of the genome retrieved and reads annotated for each threshold in sRNA abundance.
* `reads.txt` - a breakdown of the aligned reads making up the top 30% of a locus's total expression. This is useful to quickly get constituent sequences in complex loci.


### Calculation

These are secondary calculations that will be done on a tradeoff annotation.

`yasma.py count` 
This produces a file of counts for every locus. This is broken out into a separate module because yasma does this very thoroughly, producing a simple count matrix (useful for DEseq) and also a large, long format count with separates counts by locus, strand, size, and readgroup. This can save some headaches in later analyses.

`yasma.py context`
Compares loci locations with an NCBI-formatted `.gff3` file provided. Gives overlaps and nearby genes for each sRNA locus.

`yasma.py jbrowse`
This module was made to take some of the headache out of making nicely-formatted jbrowse-ready coverage and annotation maps. This produces `.bw` files for all specified sizes and strands of sRNAs, which can then be plotted in the same track with provided configuration code. If provided with a `-j, --jbrowse_directory`, this will automatically look for a config file, update it, and copy all relevant files to a directory based on the genome name.

`yasma.py hairpin` 
This tool is still in development, but it is meant to evaluate all loci for the possibility that they are derived from an RNA hairpin, rather than RDR-dsRNA. This is not finalized, but it generally looks for stranded regions and folds them, analyzing their profile based on a battery of rules from multiple publications.



### Other stuff

`utilities` includes several other functions, most of which have been used primarily for testing. Probably not relevant as of now to wider use.

</a>

<a name="pipeline"/>

## YASMA cookbook

This describes a simple pipeline making use of 5 libraries from 2 conditions. Most inputs are described in the `yasma.py inputs` call, with the pipeline calls importing these values from the `inputs.json` file.

This most simple pipeline will provide an annotation based on all libraries and a corresponding counts file.

```
mkdir srna_analysis
cd srna_analysis

## Using the following hypothetical libraries from the corresponding conditions
# lib_1.fa -> wt
# lib_2.fa -> wt
# lib_3.fa -> mut
# lib_4.fa -> mut
# lib_5.fa -> mut

## supplying all input information for the analysis
yasma.py inputs -o . \                     # making a working directory in this folder
-ul ../untrimmed_library_directory/*.fa \  # wherever these are located
-g path_to_your_genome.fa \
-c lib_1:wt lib_2:wt lib_3:mut lib_4:mut lib_5:mut  # describing library-to-condition relationships

## basic call
yasma.py adapter
yasma.py trim
yasma.py align
```

These basic steps have now trimmed and aligned the sRNAs. Now we can annotate. First, we will try a full annotation of all conditions.

```
yasma.py tradeoff    # defaults to all conditions
yasma.py count
```

This performed annotations and counts on all conditions.

Suppose we only want to annotate using a specific condition, and count all libraries against it. We'll perform a different annotation specifying the annotation conditions (can be more than one). Count automatically quantifies all logged annotation gffs. 

In this example, this might be quite important if you expect that the mutant will change or lose loci compared to wt. Suppose the mutant results in a loss of sRNAs - this might be problematic to use as the basis of a sRNA annotation.

```
## to perform the annotation with a specific condition(s)
yasma.py tradeoff -ac wt --name wt
yasma.py count
```

Of course YASMA is quite flexible, and you may give it trimmed libraries (it should actually detect this automatically).

```
## using pre-trimmed libraries
yasma.py inputs -o . \
-tl lib_1.fa lib_2.fa lib_3.fa lib_4.fa lib_5.fa \
-g path_to_your_genome.fa \
-c lib_1:hyphae lib_2:hyphae lib_3:conidia lib_4:conidia lib_5:conidia

yasma.py align
yasma.py tradeoff
yasma.py count
```

You may even supply it with a different alignment and skip the library stuff altogether. Note, this requires that the @RG flag specifies the library names that match the supplied conditions.
```
## using an alignment as input
yasma.py inputs -o . \
-a path_to_alignment.bam \ 
-c lib_1:hyphae lib_2:hyphae lib_3:conidia lib_4:conidia lib_5:conidia

yasma.py tradeoff
yasma.py count
```


</a>

<a name="Post-filtering"/>

In an effort to reduce false positives, annotated loci are filtered to remove loci (`removed_loci.gff3`) which seem likely to be noise. This is somewhat subjective, but basically it is focused on loci which are dominated by very few unique sequences. This is defined by four metrics:

### Skew
This is basically how many reads in the annotation come from the most-abundant read. Upon publication, YASMA-tradeoff filters loci which have greater than 90% (`--max_skew 0.90`) of reads are a single sequence. This is useful to remove annotations from duplicated sequences, but might also have a negative effect on highly precisely processed miRNAs.

To save highly precise sRNAs in well-studied systems, the option `--override_skew` allows this filter to be ignored for loci meeting certain sizing requirements. 

### Complexity
Complexity is focused on the number of unique reads across genomic space. This is meant to remove loci which have over-expanded and might represent a very large genomic area with only very few unique reads. This is defined by `--min_complexity`, which is the minimum number of unique reads per 1000 nt. By default this is `10 unique reads per 1000`. Low complexity loci are often caught by the skew filter also.

### Abundance
Minimum abundance is very simple: how many reads are defining a given locus. By default `--min_abundance 50`, ensures that there are at least enough reads to assess the size profile of a given locus. Loci filtered by this might be real, but YASMA-tradeoff ignores them to avoid false-confidence. If real, you must have more depth to find this locus.

### Abundance density
In addition to unique read density (complexity), YASMA-tradeoff filters loci which are generally too large for the number of containing reads. By default, this is `--min_abundance_density 100` reads per 1000 nt. This prevents expansive barely expressed loci, which can occur in edge-cases.




