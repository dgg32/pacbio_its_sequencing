import multiprocessing
import os
from threading import Semaphore
import sys
import sourmash
import configparser

import screed
import re
import itertools
import copy

rx_size = re.compile(r'size=(\d+)')



abskews = 2
k_size = 17

min_unmatched_k_mer = 30

total_account_for = 0.95
min_percentage = 0.001


kmer_dynamic = True
show_detail = False

padded = True


config_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'chimera_kmer_multiprocess.config')

if os.path.exists(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)


    abskews = float(config['DEFAULT']['abskews'])
    k_size = int(config['DEFAULT']['k_size'])

    total_account_for = float(config['DEFAULT']['total_account_for'])
    min_percentage = float(config['DEFAULT']['min_percentage'])

    kmer_dynamic = config['DEFAULT']['kmer_dynamic']

    if kmer_dynamic not in ["True", "true"]:
        kmer_dynamic = False
    else:
        kmer_dynamic = True

    show_detail = config['DEFAULT']['show_detail']

    if show_detail not in ["True", "true"]:
        show_detail = False
    else:
        show_detail = True

    padded = config['DEFAULT']['padded']
    if padded not in ["True", "true"]:
        padded = False
    else:
        padded = True



def duo (e1, e2):

    
    e2_in_e1 = e2.intersection(e1)
    
    what_left = e2.difference(e1)
    
    return ([e2_in_e1, what_left])

in_queue = multiprocessing.Manager().Queue()
out_queue = multiprocessing.Manager().Queue()
detail_queue = multiprocessing.Manager().Queue()


writeLock = Semaphore(value = 1)

#filename = "/home/sih13/Downloads/cyano8_nonchimera.fasta"
filename = sys.argv[1]

#print (filename)

seqs = {}
names = []
sizes = []

with screed.open(filename) as seqfile:
    for read in seqfile:
        #print(read.name, read.sequence)
        search_size = rx_size.search(read.name)
        
        if search_size:

            how_many_kmer = 1000

            padded_sequence = read.sequence
            if padded == True:
                padded_sequence = "A" * k_size + read.sequence + "T" * k_size

            if kmer_dynamic == True:
                how_many_kmer = len(padded_sequence) - k_size + 1

            E2 = sourmash.MinHash(n=how_many_kmer, ksize=k_size)
            E2.add_sequence(padded_sequence)
            e2 = set(E2.get_hashes())
        
            seqs[read.name] = {"size": int(search_size.group(1)), "md5s": e2}
            
            names.append(read.name)
            sizes.append(int(search_size.group(1)))

            
sorted_name = [x for _,x in sorted(zip(sizes,names))]
sorted_size = [x for x,_ in sorted(zip(sizes,names))]

#print (sorted_size[:10])

def work():
    while True:
        
        i = in_queue.get()
        
        current_investigate_name = sorted_name[i]
        current_investigate_size = sorted_size[i]
        
        #print (i)
        
        current_list = copy.deepcopy(sorted_name)
        

        e2 = seqs[current_investigate_name]["md5s"]

        #print ("current_investigate_name:", current_investigate_name)

        del (current_list[i])

        hit_list = []
        coverages = []
        #print (current_list)

        while len(current_list) > 0 and len(e2) > k_size:
            max_percentage = 0
            max_percentage_name = ""
            what_left = set()
            #print (current_list)

            for j in range(i + 1, len(current_list)):
                if seqs[current_list[j]]["size"] / current_investigate_size > abskews:

                    # E1 = sourmash.MinHash(n=(len(seqs[current_list[j]]["seq"]) - k_size + 1), ksize=k_size)


                    # E1.add_sequence(seqs[current_list[j]]["seq"])


                    #e1 = set(E1.get_hashes())

                    e1 = seqs[current_list[j]]["md5s"]


                    duo_result = duo(e1, e2)

                    overlap = len(duo_result[0]) / len(seqs[current_investigate_name]["md5s"])

                    if overlap > max_percentage and overlap >= min_percentage:
                        max_percentage = overlap
                        max_percentage_name = current_list[j]
                        what_left = duo_result[1]

            if max_percentage_name != "":


                current_list.remove(max_percentage_name)
                e2 = what_left

                hit_list.append(max_percentage_name)
                coverages.append(max_percentage)
            else:
                break

        if show_detail == True:
            detail_queue.put([current_investigate_name, hit_list, coverages])

        #print (sum(coverages))
        if len(coverages) >= 2 and sum(coverages) >= total_account_for and len(e2) <= min_unmatched_k_mer:




            out_queue.put(current_investigate_name)
            #print (to_delete_name)
            #writeLock.release()
        
        in_queue.task_done()

for i in range(5):
    t = multiprocessing.Process(target=work)
    t.daemon = True
    t.start()


for i in range(len(sorted_name)):
    #print (sorted_name[i])
    in_queue.put(i)
    
in_queue.join()

n = 0

if show_detail == True:
    while not detail_queue.empty():
        detail = detail_queue.get()
        
        print (f"Query: {detail[0]}:\n")
        for name, percentage in zip(detail[1], detail[2]):
            print (f"{name}: {percentage}")
        print ("\n\n")
        detail_queue.task_done()
detail_queue.join()

while not out_queue.empty():
    result = out_queue.get()
    print (result)
    n += 1
    out_queue.task_done()
    #print ("\n\n")

#out_queue.join()
print (f"There are {n} potential chimeras.")
#print (kmer_dynamic)
