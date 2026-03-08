import itertools 

"""
Generates a kmer : index mapping system
For any kmer size
"""
def Build_mapping_dic(Kmer_size):
    Bases = "ATCG"

    Canonical_dic = {}
    Mapping_dic = {} 
    
    # Converts a kmer into its compliment
    Compliment = str.maketrans("ATCG", "TAGC") 
    for product in sorted(list(itertools.product(Bases, repeat = Kmer_size))):
        # Join into single kmer strings  
        Kmer = "".join(product)
        # get the compliment kmer AND reverse it 
        RC_kmer = Kmer.translate(Compliment)[::-1]
        # Find the lexographicaly smallest 
        Can_kmer = min(Kmer, RC_kmer)
        # add Canonical dictionary to the dictionary
        Canonical_dic[Kmer] = Can_kmer

    Can_kmers = sorted(set(Canonical_dic.values()))

    Canonical_dic = {Can_kmer : index for index, Can_kmer in enumerate(Can_kmers)}

    # Creating a ALL_kmer mapping system 
    # Map Kmers and its reverse compliment to the same index - this is done to speed up counting kmers
    # prevents the need of converting it to the reverse compliment and then getting the index
    for Can_kmer, index in Canonical_dic.items():
        # adding Canonical kmer : index to dictionary
        Mapping_dic[Can_kmer] = index
        # Getting the reverse compliment of the kmer
        RC_kmer = Can_kmer.translate(Compliment)[::-1]
        # Adding the Reverse Compliment kmer to the dictionary 
        Mapping_dic[RC_kmer] = index

    return Mapping_dic, Canonical_dic



"""
Converts the chromsosome raw sequence, into a 1 dimensional numpy array of the kmer indices
for example, with kmer size of 4 , the following chromosome sequnce [ATCGAC] is converted into the kmer indices [4 100 69]
kmer : indexes dictionary are made  in Build_kmer_index
INVALID KMERS ARE RETURNED AS 9999

############ WARNING ################## WARNING ################# WARNING ############
- This functions ASSUMES that the bacteria has a circular chromosome
- In the future this needs to be replaced with a dynamic function which recognises circular vs linear 

"""
import numpy as np

def Chromosome_conversion(Kmer_size, Chromosome_seq, Chromosome_len, Mapping_dic):  
    
    # My plan is jsut to convert it into a very long string of 0 --> 135
    Chr_kmer_arr = np.zeros(Chromosome_len, dtype="int16")
    # Iterate through each Chromosome position
    # Artificaly Extend chromosome sequence out by the kmer size - simulates its circular nature (very minor change - will BARELY affect the kmer counts - but still)
    Chromosome_seq_ext = Chromosome_seq + Chromosome_seq[0 : Kmer_size]
    for Pos in range(len(Chromosome_seq)):
        # extract kmer, convert to canonical, convert to index (I CAN GET RID OF CONVERTING IT TO A CANONICAL - to do later)
        Kmer = Chromosome_seq_ext[Pos : Pos + Kmer_size]
        Kmer_id = Mapping_dic.get(Kmer)
        
        # append to numpy array 
        if Kmer_id is None:
            Chr_kmer_arr[Pos] = 9999
        else : Chr_kmer_arr[Pos] = Kmer_id

    return Chr_kmer_arr



## Old funciton that simply returns the kmer counts across the Chromosome kmer array 

def Calc_counts(Chr_kmer_arr, Canonical_dic, Chunk_sizes, N_chunks):
    N_can_kmers = len(Canonical_dic)
    
    Chunks_kmer_counts = np.zeros((len(Chunk_sizes), N_chunks, N_can_kmers))

    for Chunk_size_ID, Chunk_size in enumerate(Chunk_sizes):
        # Artifically extend out the chromosome to simulte its circular nature 
        Chr_kmer_arr_ext = np.append(Chr_kmer_arr, Chr_kmer_arr[:Chunk_size])

        for Chunk_ID in range(N_chunks):
            #Random Seed - using the unextended Chromosome to ensure a valid seed location is chosen
            Seed = random.randint(0, len(Chr_kmer_arr))
            # take a slice out of the Chr_kmer_arr
            Slice = Chr_kmer_arr_ext[Seed : Seed+Chunk_size]
            # Remove invalid Kmers - which are represented as 9999
            Slice_val = Slice[Slice != 9999]
            # Collapse the kmer counts 
            Counted = np.bincount(Slice_val, minlength = 136)

        
            Chunks_kmer_counts[Chunk_size_ID][Chunk_ID] = Counted

    return Chunks_kmer_counts

