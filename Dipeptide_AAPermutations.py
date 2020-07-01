import os
import sys
import itertools

amino_acids = {
    'A' : 71.03711,         #Alanine
    'R' : 156.10111,        #Arginine
    'N' : 114.04293,        #Asparagine
    'D' : 115.02694,        #Aspartic acid
    'C' : 103.00919,        #Cysteine
    'E' : 129.04259,        #Glutamic acid
    'Q' : 128.05858,        #Glutamine
    'G' : 57.02146,         #Glycine
    'H' : 137.05891,        #Histidine
    'I' : 113.08406,        #Isoleucine
    'L' : 113.08406,        #Leucine
    'K' : 128.09496,        #Lysine
    'M' : 131.04049,        #Methionine
    'F' : 147.06841,        #Phenylalanine
    'P' : 97.05276,         #Proline
    'S' : 87.03203,         #Serine
    'T' : 101.04768,        #Threonine
    'W' : 186.07931,        #Tryptophan
    'Y' : 163.06333,        #Tyrosine
    'V' :  99.06841,        #Valine
    'a' :  85.12,           #a-aminoisobutyric acid
    'b' :  99.14,           #Isovaline
    'c' :  71.09,           #b-alanine
}
valid_aminoacids = "ARNDCEQGHILKMFPSTWYVabc"
assert len(amino_acids) == len(valid_aminoacids), "Mismatch in amino acid count"

TOLERANCE = 0.05

class Sequence:
    def __init__(self, seq, cyclic, mass):
        self.sequence = seq
        self.cyclic = cyclic
        self.mass = mass
        self.counted = False

def generate_masses(output_fp, sequence_lengths, cyclic):
    sequences = []
    for length in sequence_lengths:
        permutations = itertools.product(valid_aminoacids, repeat = length)
        for p in permutations:
            seq = "".join(p)
            mass = 0
            for aa in seq:
                mass += amino_acids[aa]
            if not cyclic:
                sequence = seq 
                mass += 18.0
                cyclic_flag = 0 
            else:
                sequence = "c" + seq
                cyclic_flag = 1
            output_fp.write("%s,%d,%d,%.3f\n" % (sequence, length, cyclic_flag, mass))
            sequences.append(Sequence(sequence, cyclic, mass))
    return sequences

def write_header(output_fp):
    output_fp.write("Sequence,Length,Cyclic,Mass\n")

def write_summary(summary_fp, linear_sequences, cyclic_sequences):
    summary_fp.write("Linear sequences:\n")
    summary_fp.write("Mass,Count,Sequences\n")
    for sequence in linear_sequences:
        count = 0
        if not sequence.counted:
            sequence.counted = True
            mass = sequence.mass
            count += 1
            sequences = [sequence.sequence]
            for compare_seq in linear_sequences:
                if not compare_seq.counted and abs(compare_seq.mass - mass) < TOLERANCE:
                    compare_seq.counted = True
                    sequences.append(compare_seq.sequence)
                    count += 1
        if count > 0:
            summary_fp.write("%.2f,%d" % (mass, count))
            for s in sequences:
                summary_fp.write(",%s" % s)
            summary_fp.write("\n")
    summary_fp.write("\n")

    summary_fp.write("Cyclic sequences:\n")
    summary_fp.write("Mass,Count,Sequences\n")
    for sequence in cyclic_sequences:
        count = 0
        if not sequence.counted:
            sequence.counted = True
            mass = sequence.mass
            count += 1
            sequences = [sequence.sequence]
            for compare_seq in cyclic_sequences:
                if not compare_seq.counted and abs(compare_seq.mass - mass) < TOLERANCE:
                    compare_seq.counted = True
                    sequences.append(compare_seq.sequence)
                    count += 1
        if count > 0:
            summary_fp.write("%.2f,%d" % (mass, count))
            for s in sequences:
                summary_fp.write(",%s" % s)
            summary_fp.write("\n")
    summary_fp.write("\n")

    # Combine them all together
    overall_sequences = linear_sequences + cyclic_sequences
    # Reset the "counted" flag
    for s in overall_sequences: s.counted = False

    summary_fp.write("Combined sequences:\n")
    summary_fp.write("Mass,Count,Sequences\n")
    for sequence in overall_sequences:
        count = 0
        if not sequence.counted:
            sequence.counted = True
            mass = sequence.mass
            count += 1
            sequences = [sequence.sequence]
            for compare_seq in overall_sequences:
                if not compare_seq.counted and abs(compare_seq.mass - mass) < TOLERANCE:
                    compare_seq.counted = True
                    sequences.append(compare_seq.sequence)
                    count += 1
        if count > 0:
            summary_fp.write("%.2f,%d" % (mass, count))
            for s in sequences:
                summary_fp.write(",%s" % s)
            summary_fp.write("\n")
    summary_fp.write("\n")

def main():
    output_filename = sys.argv[1]
    output_fp = open(output_filename, "w")
    write_header(output_fp)

    linear_sequences = generate_masses(output_fp, [2, 3], cyclic =  False)
    cyclic_sequences = generate_masses(output_fp, [2], cyclic = True)
    output_fp.close()    

    summary_filename = sys.argv[2]
    summary_fp = open(summary_filename, "w")
    write_summary(summary_fp, linear_sequences, cyclic_sequences)
    summary_fp.close()

if __name__ == "__main__":
    main()
