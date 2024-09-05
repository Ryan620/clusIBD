###  clusIBD
 clusIBD is designed to  detect IBD segments  using unphased genetic data from challenging DNA samples. Two search modes are implemented, i.e. “within-search” and “across-search”. Within-search refers to inferring IBD segments for all the pairs within a database (one input file), whereas across-searching works by inferring IBD segments between two individuals from different databases (two input files). The latter is often performed when one or more individuals of interest are searched against a large database to identify potential relative. 

### Before running clusIBD
The PLINK binary ped format (i.e., bed, bim, and fam files) is accepted as inputs. So  users should convert the input data (e.g. VCF files) into PLINK binary format.
  plink --vcf example.vcf  --double-id --make-bed --chr 1-22 --geno 0.5  --mind 0.5 --out  example

### Quickstart
clusIBD -f example -o test


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


###clusIBD output
clusIBD outputs two files, which are both a human-readable text format. The .details file outputs all segments for each pair and the .summary file shows the total numbers and lengths of IBD segments for each pair.
* .details file 
  sample1	sample2	chromosome	start	end	lengths (Mb)	IBD_type

*.summary file
  sample1 sample2	total_number	total_lengths (Mb)

## License

This project is licensed under the GPL-3.0 - see the [LICENSE](LICENSE) file for details

