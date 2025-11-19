# Installation

### 1) Best option: conda

The first and best option for installation is through bioconda, which includes all dependencies.
```
conda create -n my_env
conda activate my_env
conda install bioconda:yasma

yasma -h
```

### 2) Installing from github with a native python installation

This option requires more work, but might be more flexible if you have trouble with conda or dependencies in conda. 

This simply downloads and installs the tool to your system PATH. Extra python modules are required for Yasma modules, which you must install manually. These are lazy-loaded and some may not be necessary if you don't need that analysis. Its a good idea to use a virtual environment here, which i show with `venv` in python.

```
python3 -m venv yasma
source yasma/bin/activate

git clone https://github.com/NateyJay/YASMA.git
cd YASMA

./yasma.py -h

```

More permanent installation would likely include moving this directory to a more permanent place and adding it to your path.

```
mv YASMA /usr/local/
## adding the following line to ~/.bash_profile -> export PATH="/usr/local/YASMA:$PATH
## this can be done with
nano ~/.bash_profile
# or
echo "export PATH="/usr/local/YASMA:$PATH" >> ~/.bash_profile

source ~/.bash_profile
```

### Dependencies
Yasma makes use of many tools through wrappers as well as some non-preinstalled python modules. 

Python modules/programs:
* numpy
* click
* click-option-group
* pysam
* pprintpp
* pyBigWig (**coverage**, ***jbrowse**, **tradeoff** with -bw)
* cutadapt (only needed for **trim**)
* levenshtein (**hairpin**)
* viennarna (**hairpin**)

Other programs
* bowtie1(https://bowtie-bio.sourceforge.net/index.shtml) (**align**)
* [sra-tools](https://github.com/ncbi/sra-tools) (**download**)

#### Ubuntu/Debian systems
Nearly all dependencies and tools are available through apt.

```
sudo apt install python3-numpy python3-click python3-click-option-group python3-pysam python3-pyBigWig python3-cutadapt python3-levenshtein bowtie python3-pprintpp

## viennarna is not available through apt, so we can use pip
python3 -m pip install viennarna

## sra-tools is also not in apt. Available here:  https://github.com/ncbi/sra-tools.
```

#### MacOS
Other systems like MacOS relies on pip to download packages. Note, you may need to use --break-system-packages if you don't use a virtual environment.

```
## core modules
python3 -m pip install numpy click click-option-group pysam cutadapt pyBigWig levenshtein viennarna pprintpp

## sra-tools and bowtie1 are not found in pip, but they are in homebrew
brew install brewsci/bio/bowtie
brew install sratoolkit
```
