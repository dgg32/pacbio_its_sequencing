# pacbio_its_sequencing
DSMZ PacBio 16S_ITS sequencing evaluation script collections. The scripts here are meant to be integral part of  the PacBio 16S_ITS analytic pipeline. They perform different tasks before the data can be finally analysed in a tsv format.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites


 
```
screed

sourmash

cutadapt

```

### Installing

Scripts can be run as is without installation.

## Running kmer-based chimera check

```
python chimera_kmer_multiprocessing.py input.fasta
```
### Bookending

To make sure only sequences with both forward and reverse primers are including in the analyses, a "bookend" check can be performed. First, use cutadapt to recognize sequences with primers by:

```
cutadapt -g file:'forward.fasta'  -o forwardout1  input.fasta --minimum-length 1
cutadapt -a file:'forward.fasta'  -o forwardout2  input.fasta --minimum-length 1
cutadapt -g file:'reverse.fasta'  -o reverseout1  input.fasta --minimum-length 1
cutadapt -a file:'reverse.fasta'  -o reverseout2  input.fasta --minimum-length 1
```

And then, run the script below to get the sequence fasta:

```
python bookend.py forwardout1 forwardout2 reverseout1 reverseout2 input.fasta > result_fasta
```




## Authors

* **Sixing Huang** - *Coding*
* **Heike Freese** - *Concept and testing*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to sourmash team for providing an effective way to calculate kmers and developing the concept of mashing.
* vsearch for doing the heavy-lifting - search, chimera detection and so on
* etc

