# TFBSscan

Transcription factors play an important role in the gene regulation. To search for the regulatory regions in the genome the Transcription Factor Binding Sites Scan (TFBSscan) was developed.

## Installation

1. Clone the directory
```bash
git clone https://github.com/petrokvitka/TFBScan
```
2. Switch to the directory
```bash
cd TFBSscan
```
3. Create the needed environment from the file tfbs.yaml
```bash
conda env create --file tfbs.yaml
```
4. Activate the environment
```bash
source activate tfbs
```
5. You are ready to use the TFBS scan! To learn how to use TFBS scan, just type 
```bash
python tfbsscan.py
```

## Input files

The TFBSscan requires two input files:
* file with transcription factor motifs, which can be downloaded from [JASPAR database](http://jaspar.genereg.net/) in .MEME or .PFM format;
* a file with genome in FASTA format.

To optimize the scanning, a file in the .BED format with peaks can also be provided. The peak is a so called open region or in other words the region of interest in the ATAC-seq. If the .BED file with peaks was provided, the TFBSscan run will be way faster, as we will look for the needed information only within the regions we are interested in.

## Output files

As an output the TFBSscan produces the file in .BED format with all found binding sites for each motif provided as input. If only one motif was provided as input, only one output file will be created. For three input motifs three output files in .BED format will be outputed. If the output file is empty, no binding sites for this motif were found.

## Tools within the TFBSscan

There are currently two different tools a user can choose in between when starting the TFBSscan:
* [FIMO](http://meme-suite.org/doc/fimo.html) (Find Individual Motif Occurences) as a part of the MEME Suite. There is a possibility to pass additional parameters to FIMO. TFBSscan uses by default the next parameters for FIMO call: `--text` to print the FIMO output in a simple text file and `--no-qvalue` because the calculation of q-values is not possible while using the parameter `--text`.
* [MOODS](https://github.com/jhkorhonen/MOODS) (Motif Occurrence Detection Suite).

FIMO and MOODS need the input files in different formats (.MEME for FIMO and .PFM for MOODS). The TFBSscan allows an automatic convertion of these both formats regarding on the tool the user wants to use. 

## Customizing the TFBSscan run

There are several options to customize the TFBSscan:
* `--output_directory` to set the desired output directory, by default the output files will be saved to the path ./output/;
* `--bed_file` to provide a file with peaks of interest to search within in .BED format;
* `--use [fimo/moods]` to choose the tool (by default TFBSscan uses FIMO);
* `--clean [all/nothing/cut_motifs/fimo_output/merge_output/moods_output]` set one or more option from the list to be able to chose the temporare files from the output directory after the run. By default TFBSscan deletes all temporary files and will leave only the output files with binding sites for each motif in .BED format;
* `--cores` to set how many cores the TFBSscan is allowed to use while starting the multiprocessing. By default only 2 cores are allowed;
* `--fimo "option"` pass additional options for the FIMO run within the quotation marks. For example, a background file created with [fasta-get-markov](http://meme-suite.org/doc/fasta-get-markov.html) can be passed;
* `--p_value` set the p-value for the search. By default the p-value is 1e-4 or 0.0001;
* `--resolve_overlaps` to delete the ovelaps with lower score. By default no overlaps are deleted;
* `--moods_bg 0.25 0.25 0.25 0.25` to set hte background for MOODS. By default the standard MOODS background is used which is 0.25 0.25 0.25 0.25;
* `--hide_info` to print the information about the run only to the log-file. By default all the information is printed to the terminal.

## Example data

To try how the TFBSscan works, the test data was uploaded to the folder [example](./example). The file with 3 random motifs was downloaded from the [JASPAR database](http://jaspar.genereg.net/) in .MEME format. The fasta file [mm10](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/) with the first chromosome is used to create the file in the [example_output](./example_output) folder. The TFBSscan call used for this example is:

```bash
python tfbsscan.py -g ../chr1.fa -m example/3_motifs.meme -o example_output
```

There are three output files in .BED format in the example_output folder and a [logfile](./example_output/tfbsscan.log) containing information about the run. Each of the output .BED files contains the binding sites of the corresponding motif for the first chromosome of the mm10 genome.




