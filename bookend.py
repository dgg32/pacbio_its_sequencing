import csv
import sys
import screed

forward_1_file = sys.argv[1]
forward_2_file = sys.argv[2]
reverse_1_file = sys.argv[3]
reverse_2_file = sys.argv[4]
fasta_with_length = sys.argv[5]
#bbmapskimmer.sh in=$2 out=$samout  ref=$fasta_with_length idfilter=0.9 k=8 noheader=t threads=4 ambiguous=all

ssu = set()
lsu = set()

with screed.open(forward_1_file) as seqfile:
	for read in seqfile:
		ssu.add(read.name)

with screed.open(forward_2_file) as seqfile:
	for read in seqfile:
		ssu.add(read.name)

with screed.open(reverse_1_file) as seqfile:
	for read in seqfile:
		lsu.add(read.name)

with screed.open(reverse_2_file) as seqfile:
	for read in seqfile:
		lsu.add(read.name)

with screed.open(fasta_with_length) as seqfile:
	for read in seqfile:
		seq = read.sequence

		name = read.name

		temp_name = read.name

		has_primer = False
		if name in ssu:
			#print (">" + name + ",with_ssu" + "\n" + seq)
			temp_name += ",with_ssu"
			has_primer = True
		else:
			temp_name += ",no_ssu"
		if name in lsu:
			temp_name += ",with_lsu"
			has_primer = True
		else:
			temp_name += ",no_lsu"


		print (">" + temp_name + "\n" + seq)