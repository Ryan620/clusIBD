###  clusIBD

 clusIBD is designed to  detect IBD segments  using unphased genetic data from challenging DNA samples. Two search modes are implemented, i.e. “within-search” and “across-search”. Within-search refers to inferring IBD segments for all the pairs within a database (one input file), whereas across-searching works by inferring IBD segments between two individuals from different databases (two input files). The latter is often performed when one or more individuals of interest are searched against a large database to identify potential relative. 

### Before running clusIBD

The PLINK binary ped format (i.e., bed, bim, and fam files) is accepted as inputs. So users should convert the input data (e.g. VCF files) into PLINK binary format. For example, you may use the following command line:
```
  plink --vcf example.vcf  --double-id --make-bed --chr 1-22 --geno 0.5  --mind 0.5 --out  example
```

### Quickstart
```
./bin/clusIBD -f ./example/random_100
```
* `random_100`  is the  prefix of PLINK binary ped format. Make sure that random_100.bed, random_100.bim, and random_100.fam are all  in the same dictionary.
* Because the prefix of output file is not specified, 'out' is used by default. Once finished, two files, namely 'out.IBD.details' and 'out.IBD.summary', will be generated at the present dictionary.

### clusIBD Options:

* `-f` or `---file_prefix`
	* The prefix of bed file. This is required.
* `-F` or `--file_prefix_2`
	* The prefix of another bed file.
* `-n` or `--bin_size`
	* The number of SNPs per window. The default is -1 and the number is estimated automatically.
* `-q` or `--minQ1`
	* The minimal percentage for IBD1. The default is 0.15.
* `-Q` or `--minQ2`
	* The minimal percentage for IBD3. The default is 0.05.
* `-R` or `--fpr`
	* The  false positive rate that is allowed. The default is 0.001.
* `-s` or `--size`
	* The random size for parameter estimation. The default is 20.
* `-l` or `--min_length`
	* The minimal length for an IBD segment to be considered true. The default is 5.
* `-L` or `--IBD2_length`
	* The minimal length of IBD1 for IBD2 estimation. The default is 1000.
* `-c` or `--cpu`
	*  The number of CPU cores to be used. The default is 5.
* `-p` or `--pairs_file`
	*  File containing sample pairs, each line represents a pair.
* `-o` or `--out`
	*  The prefix of output file.If it is not specified,"out" is used.

### clusIBD output

clusIBD outputs two files, which are both a human-readable text format. The .details file outputs all segments for each pair and the .summary file shows the total numbers and lengths of IBD segments for each pair.
* .details file 
```
sample1	sample2	chr	start	end	lengths (Mb)	IBD_type
NC10m	NC7f	2	53.715181	83.134022	29.418841	IBD1
NC10m	SP3m	2	44.149442	86.440963	42.291520	IBD1
NC1m	NC4m	2	217.043942	240.662813	23.618871	IBD1
NC1m	NC4m	4	61.609072	107.57351	45.964438	IBD1
```
* .summary file
```
sample1 sample2	total_number	total_lengths (Mb)
NC1m	NC4m	9	310.816892
NC1m	NE1m	1	24.9349800
NC1m	NC7f	1	36.4197279
NC1m	NC9m	11	589.610456
```

### MIT License

Copyright (c) 2024 Ryan620

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

