#!/usr/bin/env python3
# (C) 2020-2022 by Kristoffer E. Joahnsson <kristoffer.johansson@bio.ku.dk>

import sys
import os
import argparse
from numpy import exp

__version__ = 1.00

def is_aa_one_nat(sequence, additional=""):
    for aa in sequence.upper():
        if not (aa in "ACDEFGHIKLMNPQRSTVWY" or aa in additional.upper()):
            return(False)
    return(True)

def read_fasta(filename, comment_char="#;"):
    """Flexible FASTA file reader without dependencies"""
    seq_list = []
    file_id = filename.split("/")[-1].split(".")[0]
    reading_fasta = False
    with open(filename,"r") as file_handle:
        for line in file_handle:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] in comment_char:
                continue
            if words[0][0] == ">":
                if reading_fasta:
                    seq_list.append((seq_id,seq))
                seq_id = words[0][1:]
                seq = ""
                reading_fasta = True
            elif reading_fasta:
                if len(words) > 1:
                    printf("WARNING: Found FASTA line with more white space separated fields:\n%s" % (line))
                seq = seq + words[0]
            else:
                if len(words) == 1:
                    seq = words[0]
                    if not is_aa_one_nat(seq):
                        print("ERROR: Non-FASTA single-column file should have a protein sequence in first column:")
                        print(line)
                        return(None)
                    seq_id = file_id+"%05d" % (len(seq_list))
                    seq_list.append((seq_id,seq))
                elif len(words) > 1:
                    seq = words[1]
                    if not is_aa_one_nat(seq):
                        print("ERROR: Non-FASTA multi-column file should have a protein sequence in second column:")
                        print(line)
                        return(None)
                    seq_list.append((words[0],seq))
        if reading_fasta:
            seq_list.append((seq_id,seq))
    return(seq_list)

def tile_sequence(seq, tile_len, tile_offset=1, add_ct_tile=True):
    """Cut sequence into tiles of length tile_len"""
    nres = len(seq)
    tile_begin = 0
    tile_end = tile_len
    tile_list = []
    while(tile_end <= nres):
        tile_list.append(seq[tile_begin:tile_end])
        tile_begin = tile_begin + tile_offset
        tile_end = tile_begin + tile_len
    if tile_end !=  nres+1 and add_ct_tile:
        tile_list.append(seq[-tile_len:])
    return(tile_list)

# Logistic regression model with L2 regularization (lambda=0.001)
# Trained on WT data with >50 raw counts total in 4 the bins
# Peptides classified degrons if PSI<2.2 and stable if PSI>2.8
# Parameter are per amino acid counts in peptides of length 17
model = {}
model["intersect"] = -0.89102423
model["C"] =   0.37721431
model["D"] =  -0.78986558
model["E"] =  -0.65124014
model["K"] =  -0.15518666
model["R"] =  -0.02030300
model["H"] =  -0.02110156
model["N"] =  -0.32782161
model["Q"] =  -0.17676485
model["A"] =   0.10844211
model["G"] =  -0.37594135
model["S"] =  -0.09627044
model["T"] =  -0.08533912
model["V"] =   0.43746326
model["M"] =   0.31182498
model["L"] =   0.53427787
model["I"] =   0.61465146
model["F"] =   0.52882600
model["Y"] =   0.45253658
model["W"] =   0.58693535
model["P"] =  -0.25880796

def calc_degron_prob(seq_list):
    p_list = []
    for seq in seq_list:
        # This model is only valid for peptides of length 17
        assert len(seq) == 17
        lin = model["intersect"] + sum([model[aa] for aa in seq])
        responce = 1/(1+exp(-lin))
        p_list.append(responce)
    return(p_list)


if __name__ == "__main__":
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="Predict degron profile for protein sequences of length 17 or more")
    # positional arguments                                                                                                                                  
    arg_parser.add_argument("seq_input", nargs="+", metavar="SEQ",
                            help="Sequence(s) in text or FASTA format and whitespace separated. Sequences longer than 17 amino acids are tiled")
    args = arg_parser.parse_args()

    # Read all sequences into memory
    seq_list = []
    for seq_input in args.seq_input:
        if is_aa_one_nat(seq_input.upper()):
            seq_list.append(("seq%05d" % (len(seq_list)), seq_input.upper()))
        elif seq_input[-4:].lower() in [".fas",".seq"] or seq_input[-6:].lower() == ".fasta":
            if not os.path.isfile(seq_input):
                print("ERROR: Cannot find file %s" % (seq_input))
                sys.exit(2)
            seq_list.extend(read_fasta(seq_input))
        else:
            print("ERROR: Argument %s is neither a protein sequence nor a FASTA file (.fas, .seq or .fasta)" % (seq_input))
            sys.exit(2)

    # Tile sequences and calculate degron probability
    for (seq_id,seq) in seq_list:
        tile_seq = tile_sequence(seq, 17, 1)
        tile_degron_prob = calc_degron_prob(tile_seq)
        # print("# === Protein %s ===" % (seq_id))
        # for i in range(len(tile_seq)):
        #     print("%s  %.5f" % (tile_seq[i], tile_degron_prob[i]))
        resi = 9
        for i in range(len(tile_seq)):
            print("%18s  %17s  %.5f  %1s  %4d" % (seq_id, tile_seq[i], tile_degron_prob[i], tile_seq[i][8], resi))
            resi += 1
