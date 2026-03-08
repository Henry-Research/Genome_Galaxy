def get_chunk_sizes(Chromosome_len):
    Chunk_sizes = list(range(500, int(Chromosome_len / 100), 5000))
    return Chunk_sizes


import numpy as np
import random

"""
Takes N samples from the Chr_kmer_arr across all of the chunk sizes 

Returns a four dimensional matrix :: Chunk_size --> iteration --> 136 by 136 

this entire function can be made speedier - but it is okay for now 

"""

def Calc_covariance(Chr_kmer_arr, Canonical_dic, Chunk_sizes, N_chunks, PseudoAddition, Type):
    N_can_kmers = len(Canonical_dic)
        
    Chunks_kmer_cov = np.zeros((len(Chunk_sizes), N_can_kmers, N_can_kmers))


    for Chunk_size_ID, Chunk_size in enumerate(Chunk_sizes):
        # Artifically extend out the chromosome to simulte its circular nature 
        Chr_kmer_arr_ext = np.append(Chr_kmer_arr, Chr_kmer_arr[:Chunk_size])


        Temp_chunk_storage = np.zeros((N_chunks, N_can_kmers))
        for Chunk_ID in range(N_chunks):
            #Random Seed - using the unextended Chromosome to ensure a valid seed location is chosen
            Seed = random.randint(0, len(Chr_kmer_arr))
            # take a slice out of the Chr_kmer_arr
            Slice = Chr_kmer_arr_ext[Seed : Seed+Chunk_size]
            # Remove invalid Kmers - which are represented as 9999 from the Chromosome conversion function 
            Slice_val = Slice[Slice != 9999]
            # Collapse the kmer counts 
            Counted = np.bincount(Slice_val, minlength = len(Canonical_dic))

            # Add 1 to each count if PsuedoAddition = True 
            if PseudoAddition:
                Counted = Counted + 1 

            # Convert into proportions if Type == "Proportion":
            if Type == "Proportion":
                Counted = Counted / sum(Counted) 

                
            # append it 
            Temp_chunk_storage[Chunk_ID] = Counted



        # row var set to false. Each row is an observation, and each column is a variable
        Cov = np.cov(Temp_chunk_storage, rowvar = False)
        
        Chunks_kmer_cov[Chunk_size_ID] = Cov
    
    return Chunks_kmer_cov




