from itertools import product
import random
import hashlib

def Get_valid_chunk_sizes(Chromosome_length):
    Chunk_list = []
    #fixed chunk sizes 
    Chunk_sizes = [100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000]

    # Only return chunk sizes which are smaller than 1/20th of the total chromosome size 
    for size in Chunk_sizes :
        if size < Chromosome_length / 20:
            Chunk_list.append(size)
    return Chunk_list


def Build_canonical_map(kmer_size):
    bases = "ACGT"
    mapping = {}
    complement = str.maketrans("ACGT", "TGCA")
    for t in product(bases, repeat=kmer_size):
        t_str = ''.join(t)
        rc = t_str.translate(complement)[::-1]
        canonical = min(t_str, rc)
        mapping[t_str] = canonical
    return mapping

def Get_seed(Chromosome_name, Half_chunk_size, Chunk_index, Chromosome_length, Interpolator = False):
    """
    Generates a random seed (ie: random.seed(42)) using the unique String identifies for that Chromosome/Chunk_size, Chunk_id.
    """
    if Interpolator == False: 
        Chromosome_name = Chromosome_name + "False"
    seed_string = f"{Chromosome_name}_{Half_chunk_size}_{Chunk_index}"
    seed = int(hashlib.sha256(seed_string.encode()).hexdigest(), 16) % (10**8)  # Keep it in a reasonable range
    # init the random seed 
    rng = random.Random(seed)
    # find a unique (but now reproducible) location to take the chunk out of the chromosome
    Seed_location = rng.randint(int(Half_chunk_size), int(Chromosome_length - Half_chunk_size))
    return Seed_location