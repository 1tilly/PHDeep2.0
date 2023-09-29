import pyfaidx

import numpy as np

class SequenceParser():
    def load_genome_from_fasta(input_path):
        genome = pyfaidx.Fasta(input_path)
        return genome
    
    def get_fasta_ref(chr, start, end, genome):
        if  "chr" in str(chr):
            return genome[chr][int(start)-1:int(end)].seq
        else:
            return genome[f"chr{str(chr)}"][int(start)-1:int(end)].seq

    def one_hot(sequence):
        classes = 4
        encode = {"N":0, "A":1, "C":2, "G":3, "T":4, "n":0, "a":1, "c":2, "g":3, "t":4}
        num_sequence = list(map(encode.get, sequence))
        targets = np.array([num_sequence]).reshape(-1)
        oh_encoding = np.insert(np.eye(classes), 0, [0,0,0,0], axis=0) # inserting the [0,0,0,0] case for unknown bases.
        result = oh_encoding[np.array(targets).reshape(-1)]
        return result.reshape(list(targets.shape)+[classes])

    def decode_one_hot(one_hot_sequence):
        encode = {"N":-1, "A":0, "C":1, "G":2, "T":3}
        decode = {v:k for k,v in encode.items()}
        result = []
        for seq in one_hot_sequence:
            result.append(decode[np.argmax(seq) if not np.all(np.array(seq)== np.array([0,0,0,0])) else -1])
        
        return "".join(result)
