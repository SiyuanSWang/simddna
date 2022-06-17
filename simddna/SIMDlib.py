from itertools import combinations, product
import numpy as np
import math
from Bio import pairwise2
import gzip
import copy

n_bits = 4

def all_n_subs(word, i):
    # Produces a generator
    # Be sure to list() and set() to remove duplicates
    # Returns all possible mismatches up to i subs from given seq
    # From stackoverflow #11679855
    for d in range(1, i+1):
        for locs in combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGT" if l != origChar]
            for poss in product(*thisWord):
                yield "".join(poss)   

def all_n_dels(seq, n):
    # Produces a generator
    # Be sure to list() and set() to remove duplicates.
    # All possible deletions up to n subs from given seq
    for d in range(1, n+1):
        for locs in combinations(range(len(seq)), d):
            deletion = [seq[i] for i in range(len(seq)) if i not in locs]
            yield "".join(deletion)
                
def all_n_ins(seq, n):
    # Returns set of all sequences with up to nins insertions from given seq.
    # From github/hawkjo/champ/seqtools.py
    outset = set()
    bases = "ACGT"
    for nins in range(1, n+1):
        for tup in combinations(range(1, len(seq) + 1), r=nins):
            for ins_bases in product(bases, repeat=nins):
                assert len(ins_bases) == len(tup), (tup, ins_bases)
                newseq = seq[:tup[0]]
                for base_idx, (i, j) in enumerate(zip(tup, tup[1:])):
                    newseq += ins_bases[base_idx] + seq[i:j]
                newseq += ins_bases[-1] + seq[tup[-1]:]
                assert len(newseq) == len(seq) + nins, (tup, newseq)
                outset.add(newseq)
    return outset

def dna_complement(N):
    # N can be a single nucleotide or sequence
    # Gives complement back in 5' -> 3'

    # INPUT: sequence (string)
    # OUTPUT: complement (string)
    
    comps = {'A':'T','T':'A','C':'G','G':'C'}
    complement = ''
    for n in N:
        complement += comps[n]
    return(complement[::-1])

def most_common(List):
    count = {}
    for item in List:
        if item in count.keys():
            count[item] += 1
        else: 
            count[item] = 1
    
    word_count = []
    for key, value in count.items():
        word_count.append((key, value))
    word_count = sorted(word_count, key=lambda x: x[1])
    
    max_occur = max([x[1] for x in word_count])
    most_common = []
    for key, value in word_count:
        if value == max_occur: 
            most_common.append(key)
    
    return(most_common)

def barcode_identifier(ReadPairDict, OrderedFwdBarcodes, OrderedRevBarcodes, \
                       barcode_distance=238):
    # Identifies the forward and reverse barcodes by searching through the whole sequence
    # Probably an inefficient way to check for where the barcode is, but more straightforward
    
    # Set the bounds for where to look for 2nd barcode
    barcode_length = len(OrderedFwdBarcodes[0][0])
    start = barcode_distance + math.floor(barcode_length/2)
    ### Note that this might be to blame for low barcode calls!
    end = start + 2*barcode_length ###

    
    fwd_key, rev_key = sorted(ReadPairDict.keys())
    fwd_seq = ReadPairDict[fwd_key][0]
    rev_seq = ReadPairDict[rev_key][0]
    
    fwdID = {"F":[], "R":[]} # which read showed the barcode 
    fwdID_idx = {"F":[], "R":[]} # position counting from the 5' end of the read that has barcode
    revID = {"F":[], "R":[]} 
    revID_idx = {"F":[], "R":[]} 
    
    for direction, seq, all_barcodes in zip(["F","R"],[fwd_seq, rev_seq], \
                                            [OrderedFwdBarcodes, OrderedRevBarcodes]):
        seq_1 = seq[:int(len(seq)/2)]
        seq_2 = seq[int(len(seq)/2):]

        for i, barcode_set in enumerate(all_barcodes):
            for barcode in barcode_set:
                if barcode in seq_1:
                    if direction == "F": # flip if reverse read, bc first barcode would be reverse
                        # locate what index it was found
                        fwdID_idx[direction].append(seq.index(barcode))
                        fwdID[direction].append(i+1)
                    else:
                        revID_idx[direction].append(seq.index(barcode))
                        revID[direction].append(i+1)
                    break
        if direction == "F":
            if len(fwdID[direction]) == 0: # Didn't find it
                fwdID_idx[direction].append(None)
                fwdID[direction].append(0)
        else: # reverse read; didn't find fwd barcode
            if len(revID[direction]) == 0: # At the last barcode, still haven't found it
                revID_idx[direction].append(None)
                revID[direction].append(0)
                    
        # Second half, look for reverse barcodes, narrow down a range
        if direction == "F":
            if fwdID[direction][0] != 0: # Fwd barcode was found; narrow range to find reverse
                # second seq is a list of sequences to look through, ordered by former barcode
                seq_2_list = []
                for idx in fwdID_idx[direction]:
                    seq_2_list.append(seq[start+idx:end+idx])
            else: # didn't find first barcode
                seq_2_list = [seq_2]
        else: # reverse
            if revID[direction][0] != 0: # Rev barcode was found; narrow range to find fwd
                # second seq is a list of sequences to look through, ordered by former barcode
                seq_2_list = []
                for idx in revID_idx[direction]:
                    seq_2_list.append(seq[start+idx:end+idx])
            else: # didn't find first barcode
                seq_2_list = [seq_2]

        for k, segment in enumerate(seq_2_list):
            for i, barcode_set in enumerate(all_barcodes[:4]):
                for j, barcode in enumerate(barcode_set):
                    if barcode in segment: 
                        if direction == "F":
                            revID_idx[direction].append(seq_2.index(barcode) + len(seq_1))
                            revID[direction].append(i+1)
                        else:
                            fwdID_idx[direction].append(seq_2.index(barcode) + len(seq_1))
                            fwdID[direction].append(i+1)
                        break
            if direction == "F":
                if len(revID[direction]) < k+1: # At the last barcode, still haven't found it
                    revID_idx[direction].append(None)
                    revID[direction].append(0)
            else: # reverse read; didn't find fwd barcode
                if len(fwdID[direction]) < k+1: # At the last barcode, still haven't found it
                    fwdID_idx[direction].append(None)
                    fwdID[direction].append(0)
    
    # Barcode-calling - determine what we're going to call the fwd and rev barcodes
    # If one is empty, default to the other
    # See if a common barcode exists; also the most commonly cited barcode
    # If there's a disagreement or neither determined it, call it an unknown barcode
    
    final_barcode_call = []
    for barcode_ID, barcode_idx in zip([fwdID, revID], [fwdID_idx, revID_idx]):
        candidates = barcode_ID["F"] + barcode_ID["R"]
        candidates = [x for x in candidates if x != 0]
        if len(candidates) != 0:
            possible_bc = most_common(candidates)
            if len(possible_bc) > 1: # conflicting barcode calls; can't report
                final_barcode_call.append(None)
            else: # choices narrowed to one
                final_barcode_call.append(possible_bc[0])
        else: # all went to 0
            final_barcode_call.append(None)
    
    return(final_barcode_call)

def n_noncontig_gaps(seq, n_gaps):
    # Given a sequence, gives a list of all possible placements of n_gaps 
    # such that they're a) not contiguous b) not on the same side
    # If n_gaps is even, will distribute gaps to either half
    
    gapped_seqs = []
    
    if n_gaps > 1: # Pair up all the gaps
        if n_gaps % 2 == 0: # even number; split between first and second half
            first_half = seq[:math.ceil(len(seq)/2)]
            second_half = seq[math.floor(len(seq)/2):]
            first_gaps = combinations(np.arange(0,len(first_half)), int(n_gaps/2))
            second_gaps = combinations(np.arange(len(first_half), len(seq)+1), int(n_gaps/2))
            possible_gaps = [x + y for x, y in product(first_gaps, second_gaps)]
            
        else: # Distribute normally
            possible_gaps = combinations(np.arange(0, len(seq)+1), n_gaps)
            
        for gaps in possible_gaps:
            gap_pos_pairs = []
            gaps = sorted(gaps) # might be unnecessary?
            for i, gap in enumerate(gaps[:-1]):
                gap_pos_pairs.append((gap, gaps[i+1]))
            # Combine unmutated segments
            inner_gaps = "-".join([seq[x:y] for x,y in gap_pos_pairs])
            all_gaps = seq[:gap_pos_pairs[0][0]] + "-" + inner_gaps + "-" + seq[gap_pos_pairs[-1][-1]:]
            gapped_seqs.append(all_gaps)
                
    elif n_gaps == 1: # Only one location to put gap
        for gaps in combinations(np.arange(0, len(seq)+1), n_gaps):
            gapped_seqs.append(seq[:gaps[0]] + "-" + seq[gaps[0]:])

    else: # n_gaps = 0
        gapped_seqs.append(seq)
            
    return(gapped_seqs)

def cell_bit_caller(fwd_seq, rev_seq, fwd_cells, rev_cells, fwd_bits_gapped, rev_bits_gapped, \
            algorithm="counting", read_direction='both', verbose=False):
    # Finds the value of the bits for each cell using results from pairwise2 alignment
    # fwd and rev cells are the unmutated constant regions for each cell
    # fwd and rev bits_gapped are the gapped variations on the 11-nt region with the bit nucleotide
    # direction 
    # algorithm changes the bit values, to whether mismatch = 1 or 0
    
    if algorithm == "counting":
        matchValue = "0"
        mismatchValue = "1"
    elif algorithm == "rule110":
        matchValue = "1"
        mismatchValue = "0"
    
    # Align fwd and rev using larger regions
    bit_id = []
    for i, (fwd_cell, rev_cell) in enumerate(zip(fwd_cells, rev_cells)):
        if read_direction in ("forward", "both"):
            fwd_align = pairwise2.align.localms(fwd_cell, fwd_seq, 1, -0.5, -0.5, -0.5)
            _, _, fwd_score, _, _ = fwd_align[0]
            fwd_result_str = str(pairwise2.format_alignment(*fwd_align[0]))
            fwd_aligned_ref = fwd_result_str.split("\n")[0].split()[1]
            fwd_aligned_seq = fwd_result_str.split("\n")[2].split()[1]
            
        if read_direction in ("reverse", "both"):
            rev_align = pairwise2.align.localms(rev_cell, rev_seq, 1, -0.5, -0.5, -0.5)
            _, _, rev_score, _, _ = rev_align[0]
            rev_result_str = str(pairwise2.format_alignment(*rev_align[0]))
            rev_aligned_ref = rev_result_str.split("\n")[0].split()[1]
            rev_aligned_seq = rev_result_str.split("\n")[2].split()[1]
        
        # Compare the two scores, use it to narrow down where to look for bit values
        if read_direction == "both":
            if fwd_score >= rev_score:
                # If they're the same, just pick fwd
                ref_seq = fwd_aligned_ref
                aligned_seq = fwd_aligned_seq
                selected_bits = fwd_bits_gapped[i]
            else:
                ref_seq = rev_aligned_ref
                aligned_seq = rev_aligned_seq
                selected_bits = rev_bits_gapped[i]

        elif read_direction == "forward": 
            ref_seq = fwd_aligned_ref
            aligned_seq = fwd_aligned_seq
            selected_bits = fwd_bits_gapped[i]

        elif read_direction == "reverse": # Use the reverse
            ref_seq = rev_aligned_ref
            aligned_seq = rev_aligned_seq
            selected_bits = rev_bits_gapped[i]
        
        # Use list of gapped bit seqs to find bit values
        for bit in selected_bits:
            # Look for any instance of this bit region
            if bit in ref_seq: 
                start = ref_seq.index(bit)
                sample_bit_seq = aligned_seq[start: start+len(bit)]
                sample_bit_seq = "".join(sample_bit_seq.split("-"))
                sample_bit_nt = sample_bit_seq[math.floor(len(sample_bit_seq)/2)]
                ref_bit_seq = ref_seq[start: start+len(bit)]
                ref_bit_seq = "".join(ref_bit_seq.split("-"))
                ref_bit_nt = ref_bit_seq[math.floor(len(ref_bit_seq)/2)]
                
                # Check if they're the same
                if sample_bit_nt == ref_bit_nt:
                    bit_id.append(matchValue)
                else: # Bits don't agree
                    # All mismatch values ought to change to G for Fwd and C for Rev
                    if read_direction == "both":
                        if (fwd_score >= rev_score and sample_bit_nt == "G") or \
                            (fwd_score < rev_score and sample_bit_nt == "C"): # Rev
                            bit_id.append(mismatchValue)
                        else:
                            bit_id.append("?")
                    elif read_direction == "forward" and sample_bit_nt == "G":
                        bit_id.append(mismatchValue)
                    elif read_direction == "reverse" and sample_bit_nt == "C":
                        bit_id.append(mismatchValue)
                    else: # Don't have the correct mutation; it's invalid
                        bit_id.append("?")
                
                break # Already found bit, shouldn't occur twice
        
        if len(bit_id) == i and verbose: # Didn't find a match to bit, this region probably didn't qualify the read
            print("Didn't find bit region for cell {}".format(i+1))
        if len(bit_id) == i: # Actually patch the hole
            bit_id.append("?")
        
    return(bit_id)

def cell_bit_caller_phred(seqs_dict, cells, bits_gapped, \
            algorithm="counting", max_mm=5, verbose=False):
    # New version of cell bit caller that uses Phred score and more stringent
    # criteria to consider reads
    # Finds the value of the bits for each cell using results from pairwise2 alignment
    # fwd and rev cells are the unmutated constant regions for each cell
    # fwd and rev bits_gapped are the gapped variations on the 11-nt region with the bit nucleotide
    # direction 
    # algorithm changes the bit values, to whether mismatch = 1 or 0
    
    if algorithm == "counting":
        matchValue = "0"
        mismatchValue = "1"
    elif algorithm == "rule110":
        matchValue = "1"
        mismatchValue = "0"

    fwd_key, rev_key = sorted(seqs_dict.keys())
    fwd_seq, fwd_Phred = seqs_dict[fwd_key]
    rev_seq, rev_Phred = seqs_dict[rev_key]
    fwd_cells, rev_cells = cells
    fwd_bits_gapped, rev_bits_gapped = bits_gapped
    
    # Align fwd and rev using larger regions
    bit_id = {"F":[],"R":[]} 
    aligned_idx = {"F":[],"R":[]}
    for i, (fwd_cell, rev_cell) in enumerate(zip(fwd_cells, rev_cells)):
        fwd_align = pairwise2.align.localms(fwd_cell, fwd_seq, 1, -0.5, -0.5, -0.5)
        _, _, fwd_score, _, _ = fwd_align[0] 
        fwd_result_str = str(pairwise2.format_alignment(*fwd_align[0]))
        fwd_aligned_ref = fwd_result_str.split("\n")[0].split()[1]
        fwd_aligned_seq = fwd_result_str.split("\n")[2].split()[1]
        fwd_aligned_idx = fwd_result_str.split("\n")[2].split()[0]
            
        rev_align = pairwise2.align.localms(rev_cell, rev_seq, 1, -0.5, -0.5, -0.5)
        _, _, rev_score, _, _ = rev_align[0] 
        rev_result_str = str(pairwise2.format_alignment(*rev_align[0]))
        rev_aligned_ref = rev_result_str.split("\n")[0].split()[1]
        rev_aligned_seq = rev_result_str.split("\n")[2].split()[1]
        rev_aligned_idx = rev_result_str.split("\n")[2].split()[0]
        
        fwd_selected_bits = fwd_bits_gapped[i]
        rev_selected_bits = rev_bits_gapped[i]

        # Use list of gapped bit seqs to find bit values
        for direction, selected_bits, ref_seq, aligned_seq, idx, align_score in \
                zip(["F","R"],[fwd_selected_bits,rev_selected_bits], \
                    [fwd_aligned_ref, rev_aligned_ref], [fwd_aligned_seq,rev_aligned_seq], \
                    [fwd_aligned_idx, rev_aligned_idx], [fwd_score, rev_score]):
            
            # Alignment must have at most max_mm dashes, or mismatch/gaps
            # Don't go by align score bc there can be perfectly matched subsets
            if aligned_seq.count("-") > max_mm or ref_seq.count("-") > max_mm:
                bit_id[direction].append("?")
                aligned_idx[direction].append(None)
                continue

            for bit in selected_bits:
                # Look for any instance of this bit region
                if bit in ref_seq: 
                    start = ref_seq.index(bit)
                    sample_bit_seq = aligned_seq[start: start+len(bit)]
                    sample_bit_seq = "".join(sample_bit_seq.split("-"))
                    sample_bit_nt = sample_bit_seq[math.floor(len(sample_bit_seq)/2)] # middle nt
                    ref_bit_seq = ref_seq[start: start+len(bit)]
                    ref_bit_seq = "".join(ref_bit_seq.split("-"))
                    ref_bit_nt = ref_bit_seq[math.floor(len(ref_bit_seq)/2)]
                    
                    phredIndex = int(idx) + start + math.floor(len(ref_bit_seq)/2)-\
                                aligned_seq[:start].count('-')-1

                    # Check if they're the same
                    if sample_bit_nt == ref_bit_nt:
                        bit_id[direction].append(matchValue)
                        aligned_idx[direction].append(phredIndex)
                    else: # Bits don't agree
                        # All mismatch values ought to change to G for Fwd and C for Rev
                        if (direction == "F" and sample_bit_nt == "G") or \
                            (direction == "R" and sample_bit_nt == "C"): # Rev
                            bit_id[direction].append(mismatchValue)
                            aligned_idx[direction].append(phredIndex)
                        else:
                            bit_id[direction].append("?")
                            aligned_idx[direction].append(None)

                    break # Already found bit, shouldn't occur twice
        
            # If not available, fill in ?
            if len(bit_id[direction]) == i : 
                # Didn't find a match to bit, this region probably didn't qualify the read
                bit_id[direction].append("?")
                aligned_idx[direction].append(None)
                if verbose:
                    print("Didn't find bit region for cell {} in read {}".format(i+1, direction))
                
    bitcall = []
    # Pick which value to go with
    for i in range(len(fwd_cells)):
        if bit_id["F"][i] == bit_id["R"][i]: # Both reads agree
            bitcall.append(bit_id["F"][i])
        # If not, print something, and go with the higher phred score
        else:
            if aligned_idx["R"][i] == None:
                bitcall.append(bit_id["F"][i])
            elif aligned_idx["F"][i] == None:
                bitcall.append(bit_id["R"][i])
            elif fwd_Phred[aligned_idx["F"][i]] >= rev_Phred[aligned_idx["R"][i]]:
                bitcall.append(bit_id["F"][i])
            else:
                bitcall.append(bit_id["R"][i])
        
    return(bitcall)


def ngs_fastq_reader(filename):
# Generator that takes an ngs fastq file and gives you data read by read
    with gzip.open(filename, "r") as f:
        file_content = f.read()
        file_content = file_content.decode("utf-8")
        file_by_line = file_content.splitlines()
        
        # Go through the lines
        for readstart_index in range(int(len(file_by_line)/4)):
            # Every 4 lines is data on a single read
            readname = file_by_line[readstart_index*4] 
            read_id, PE_id = readname.split(" ")
            seq = file_by_line[readstart_index*4+1]
            # +
            Qscore = file_by_line[readstart_index*4+3]
            
            yield (read_id, PE_id, seq, Qscore)


# A metric for looking at correlated bits - expected error for each initial value
def hamming_dist(seq1, seq2):
    # Check that they're the same length
    if len(seq1) != len(seq1):
        print("Sequences are different lengths!")
        return
    
    score = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            score += 1
    return(score)

def barcode_identifier_JA20290(ReadPairDict, OrderedFwdBarcodes, OrderedRevBarcodes, \
                       pred_second_range_fwd = (190, 231), pred_second_range_rev = (195, 260), \
                       barcode_distance=238, read_direction='both', unanimous=False):
    # Identifies the forward and reverse barcodes by searching through the whole sequence
    # V2: takes best estimate of where the barcode is located in the Reverse and Forward reads
    # pred_second_range_fwd: min and max of range, where we look for the second barcode of fwd read
    # same for pred_second_range_rev, but for reverse read
    # Index is relative to 5' -> 3' directional indexing of reads
    # For JA20290, the reverse reads only catch 1 barcode
    # New optional argument: read_direction (which reads to consider, can be both, forward, or reverse)
    # unanimous option: only makes barcode call if that barcode is the only candidate
    
    # Set the bounds for where to look for 2nd barcode
    barcode_length = len(OrderedFwdBarcodes[0][0])
    
    fwd_key, rev_key = sorted(ReadPairDict.keys())
    fwd_seq = ReadPairDict[fwd_key][0]
    rev_seq = ReadPairDict[rev_key][0]
    
    fwdID = {"F":[], "R":[]} # which read showed the barcode 
    fwdID_idx = {"F":[], "R":[]} # position counting from the 5' end of the read that has barcode
    revID = {"F":[], "R":[]} 
    revID_idx = {"F":[], "R":[]} 
    
    # Forward and Reverse are treated differently now
    for direction, seq, all_barcodes in zip(["F", "R"],[fwd_seq, rev_seq], \
                                            [OrderedFwdBarcodes, OrderedRevBarcodes]):
        if direction == "F":
            seq_1 = seq[0:2*barcode_length] # This is now narrowed down to a specific region
        else:
            seq_1 = seq[pred_second_range_rev[0]: pred_second_range_rev[1]] 
        for i, barcode_set in enumerate(all_barcodes):
            for barcode in barcode_set: # All degenerate barcodes for one OG barcode
                if barcode in seq_1:
                    # locate what index it was found
                    if direction == "F": fwdID_idx[direction].append(seq_1.index(barcode))
                    else: fwdID_idx[direction].append(seq_1.index(barcode)+pred_second_range_rev[0])
                    
                    fwdID[direction].append(i+1)
                    break
                    
        # Only for the forward reads, attempt to find a second barcode
        if direction == "F":
            # Determine a range to look for possible 2nd barcodes
            if len(fwdID[direction]) != 0: # Have found fwd barcodes
                seq_2_range = [1000, 0] # Min max range to look at
                for idx in fwdID_idx[direction]:
                    lower_idx, upper_idx = idx, idx+barcode_length
                    second_barcode_lower_idx = barcode_distance + lower_idx + math.floor(barcode_length/2)
                    seq_2_range[0] = min(seq_2_range[0], second_barcode_lower_idx)
                    seq_2_range[1] = max(seq_2_range[1], second_barcode_lower_idx + 2*barcode_length)
            else: # didn't find first barcode
                seq_2_range = list(pred_second_range_fwd)
                
            seq_2 = seq[seq_2_range[0]:seq_2_range[1]]
            for i, barcode_set in enumerate(all_barcodes):
                for j, barcode in enumerate(barcode_set):
                    if barcode in seq_2: 
                        revID_idx[direction].append(seq_2.index(barcode) + seq_2_range[0])
                        revID[direction].append(i+1)
                        break
                        
    # Barcode-calling - determine what we're going to call the fwd and rev barcodes
    # If one is empty, default to the other
    # See if a common barcode exists; also the most commonly cited barcode
    # If there's a disagreement or neither determined it, call it an unknown barcode
    
    final_barcode_call = []
    for barcode_ID, barcode_idx in zip([fwdID, revID], [fwdID_idx, revID_idx]):
        if read_direction == "both":
            candidates = barcode_ID["F"] + barcode_ID["R"] 
            candidates = [x for x in candidates]

        elif read_direction == "forward":
            candidates = barcode_ID["F"]
            candidates = [x for x in candidates]

        elif read_direction == "reverse": 
            candidates = barcode_ID["R"]
            candidates = [x for x in candidates]
        
        if len(candidates) != 0:
            if unanimous:
                candidate = candidates[0]
                added = False
                for vote in candidates:
                    if vote != candidate:
                        final_barcode_call.append(None)
                        added = True
                        break
                if not added: final_barcode_call.append(candidate)
            else:    
                possible_bc = most_common(candidates)
                if len(possible_bc) > 1: # conflicting barcode calls; can't report
                    final_barcode_call.append(None)
                else: # choices narrowed to one
                    final_barcode_call.append(possible_bc[0])
        else: # all went to 0
            final_barcode_call.append(None)
    
    return(final_barcode_call)

def csv_to_sanger_data(filename):
    loaded_data =  np.loadtxt("sanger_quantification/{}".format(filename), dtype="str",
                              delimiter=",", skiprows=2, usecols=range(2,10))
    fwd_data = []
    rev_data = []
    for j, row in enumerate(loaded_data):
        n_data = []
        for i in range(4):
            bit_0, bit_1 = row[2*i], row[2*i+1]
            try:
                bit_0, bit_1 = int(bit_0), int(bit_1)
            except ValueError:
                pass
            n_data.append((bit_0, bit_1))

        if j % 2 == 0: # odd, bc 0-index
            fwd_data.append(n_data)
        else:
            rev_data.append(n_data)
            
    return(fwd_data, rev_data)

def NGS_barcode_identifier_RAMSIMD(readPairDict, orderedFwdBarcodes, orderedRevBarcodes,\
                       pred1stBarcodeRange=(0,20), pred2ndBarcodeRange=(230,263),\
                        barcodeDistance=231, barcodeLength=11):
    # Updated to suit new design
    # Given a single pair of readnames/reads, identifies the front and back barcodes
    # These reads look to be symmetrical, so the ranges should work in both directions
    
    fwdKey, revKey = sorted(readPairDict.keys())
    fwdSeq = readPairDict[fwdKey][0]
    revSeq = readPairDict[revKey][0]
    
    fiveID = {"F":[], "R":[]} # which read showed the barcode 
    fiveIDIndex = {"F":[], "R":[]} # position counting from the 5' end of the read that has barcode
    threeID = {"F":[], "R":[]} 
    threeIDIndex = {"F":[], "R":[]} 
    
    for direction, seq, allBarcodes in zip(["F", "R"], [fwdSeq, revSeq], \
                                [orderedFwdBarcodes, orderedRevBarcodes]):
        
        # Find the 5' (aka "forward") NGS barcode
        barcode1Region = seq[pred1stBarcodeRange[0]:pred1stBarcodeRange[1]]
        # Code is much simpler since we can expect reads to be symmetrical    
        for i, barcodeSet in enumerate(allBarcodes):
            for barcode in barcodeSet:
                if barcode in barcode1Region:
                    if direction == "F": # Got this on forward read, so this is 5' barcode
                        fiveIDIndex[direction].append(barcode1Region.index(barcode))
                        fiveID[direction].append(i+1)
                    else: # First barcode from the reverse read, so it's the 3' barcode
                        threeIDIndex[direction].append(barcode1Region.index(barcode))
                        threeID[direction].append(i+1)
                    break # Found this barcode; see if there are more
                
        # Now try to find second barcode (or first, if didn't find that)
        barcode2LowerBound, barcode2UpperBound = pred2ndBarcodeRange
        
        # For forward reads
        if direction == "F": 
            if len(fiveID[direction]) != 0:
                # Found the first barcode, use it to determine where to look next
                for index in fiveIDIndex[direction]: # update region if 5' found
                    newLowerBound = index + barcodeDistance - barcodeLength
                    newUpperBound = index + 2*barcodeLength + barcodeDistance
                    if newLowerBound < barcode2LowerBound:
                        barcode2LowerBound = newLowerBound
                    if newUpperBound > barcode2UpperBound:
                        barcode2UpperBound = newUpperBound
            barcode2Region = seq[barcode2LowerBound:barcode2UpperBound]
            
            # Now try to identify it
            for i, barcodeSet in enumerate(allBarcodes):
                for barcode in barcodeSet:
                    if barcode in barcode2Region:
                        threeIDIndex[direction].append(barcode2Region.index(barcode) + \
                                                       barcode2LowerBound)
                        threeID[direction].append(i+1)
                        break
        
        # Repeat for reverse reads, which find 3' barcode first
        else: 
            if len(threeID[direction]) != 0:
                for index in threeIDIndex[direction]:
                    newLowerBound = index + barcodeDistance - barcodeLength
                    newUpperBound = index + 2*barcodeLength + barcodeDistance
                    if newLowerBound < barcode2LowerBound:
                        barcode2LowerBound = newLowerBound
                    if newUpperBound > barcode2UpperBound:
                        barcode2UpperBound = newUpperBound
            barcode2Region = seq[barcode2LowerBound:barcode2UpperBound]
            
            # Now try to identify it
            for i, barcodeSet in enumerate(allBarcodes):
                for barcode in barcodeSet:
                    if barcode in barcode2Region:
                        fiveIDIndex[direction].append(barcode2Region.index(barcode) + \
                                                       barcode2LowerBound)
                        fiveID[direction].append(i+1)
                        break
            
    # Actually call what the barcodes were
    finalBarcodeCall = []
    for barcodeID in [fiveID, threeID]:
        candidates = barcodeID["F"] + barcodeID["R"]
        ###candidates = [x for x in candidates]
        if len(candidates) != 0:
            possibleBarcode = most_common(candidates)
            if len(possibleBarcode) > 1: # more than 1 barcode is the most common
                finalBarcodeCall.append(None)
            else: # only one
                finalBarcodeCall.append(possibleBarcode[0])
        else: # No barcodes identified
            finalBarcodeCall.append(None)

    return(finalBarcodeCall)

def SIMD_barcode_identifier_full_coverage(readPairDict, orderedFwdBarcodes, orderedRevBarcodes,\
                       predRangeinFwdRead=(181,241), predRangeinRevRead=(11,71)):
    
    fwdKey, revKey = sorted(readPairDict.keys())
    fwdSeq = readPairDict[fwdKey][0][predRangeinFwdRead[0]:predRangeinFwdRead[1]]
    revSeq = readPairDict[revKey][0][predRangeinRevRead[0]:predRangeinRevRead[1]]
    
    # What the two reads think the barcode should be
    barcodeFwdRead = 0
    barcodeRevRead = 0
    
    # Get the barcode from forward read
    for i, barcodeSet in enumerate(orderedFwdBarcodes):
        for barcode in barcodeSet:
            if barcode in fwdSeq:
                barcodeFwdRead=i+1
                break
        if barcode in fwdSeq: break
    
    # Get the barcode from reverse read
    for i, barcodeSet in enumerate(orderedRevBarcodes):
        for barcode in barcodeSet:
            if barcode in revSeq:
                barcodeRevRead=i+1
                break
        if barcode in revSeq: break
        
    return(barcodeFwdRead, barcodeRevRead)

def SIMD_barcode_identifier_incomplete_coverage(readPairDict, \
                orderedFwdBarcodes, orderedRevBarcodes,\
                predRangeinRevRead=(21,46)): 
    # SIMD barcode is only expected to be in the reverse read
    # in region 26-41, +5nts buffer
    # Original settings are for Reg8
    
    fwdKey, revKey = sorted(readPairDict.keys())
    revSeq = readPairDict[revKey][0][predRangeinRevRead[0]:predRangeinRevRead[1]]
    
    # Slightly different approach from before: see if multiple are found
    barcodeRevRead = []
    
    # Get the barcode from reverse read
    for i, barcodeSet in enumerate(orderedRevBarcodes):
        for barcode in barcodeSet:
            if barcode in revSeq:
                barcodeRevRead.append(i+1)
                break
    
    if len(barcodeRevRead) == 0: # No barcode identified
        return(0)
    
    finalBarcodeCall = most_common(barcodeRevRead)
    if len(finalBarcodeCall) > 1:
        print("More than one SIMD barcode identified: ", finalBarcodeCall)
        return(0)
    else: 
        return(finalBarcodeCall[0])


def csv_to_barcode_data(filename):
    loaded_data =  np.loadtxt("sanger_quantification/{}".format(filename), dtype="str",
                              delimiter=",", skiprows=2, usecols=10)
    fwd_barcodes = []
    rev_barcodes = []
    for i, barcode in enumerate(loaded_data):
        if i%2 == 0: # Fwd
            fwd_barcodes.append(barcode)
        else: # Rev
            rev_barcodes.append(barcode)
    
    return(fwd_barcodes, rev_barcodes)

def NGS_barcode_identifier_incomplete_coverage(readPairDict, \
            orderedFwdBarcodes, orderedRevBarcodes,\
            pred1stBarcodeRange=(0,13), barcodeLength=11):
    # Version of the NGS barcode identifier for products that aren't fully
    # read, i.e. the read length is shorter than the target region
    # We only expect 1 barcode to be in each read
    
    fwdKey, revKey = sorted(readPairDict.keys())
    fwdSeq = readPairDict[fwdKey][0]
    revSeq = readPairDict[revKey][0]
    
    # For five prime/FWD barcode or three prime/REV barcode
    fiveID = [] 
    threeID = []
    
    for direction, seq, allBarcodes in zip(["F", "R"], [fwdSeq, revSeq], \
                                [orderedFwdBarcodes, orderedRevBarcodes]):
        
        # Find the barcode towards the 5' end
        barcodeFwdRegion = seq[pred1stBarcodeRange[0]:pred1stBarcodeRange[1]]
        # Code is much simpler since we can expect reads to be symmetrical    
        for i, barcodeSet in enumerate(allBarcodes):
            for barcode in barcodeSet:
                if barcode in barcodeFwdRegion:
                    if direction == "F": 
                        fiveID.append(i+1)
                    else: 
                        threeID.append(i+1)
                    break 
                
        # No need to try and find the reverse barcode; only expect 1 barcode
            
    # Actually call what the barcodes were
    finalBarcodeCall = []
    for candidates in [fiveID, threeID]:
        if len(candidates) != 0:
            possibleBarcode = most_common(candidates)
            if len(possibleBarcode) > 1: # more than 1 barcode is the most common
                finalBarcodeCall.append(None)
            else: # only one
                finalBarcodeCall.append(possibleBarcode[0])
        else: # No barcodes identified
            finalBarcodeCall.append(None)

    return(finalBarcodeCall)

def NGS_barcode_identifier_full_coverage(readPairDict, \
            orderedFwdBarcodes, orderedRevBarcodes,\
            pred5BarcodeRange=(0,13), pred3BarcodeRange=(228,249),\
            barcodeLength=11):
    # New default version of the NGS barcode identifier 
    # for products that are fully read, i.e. read length 
    # covers entire target region
    # We expect both barcodes to be in each read
    # Assumes read symmetry for barcode position!
    # Assumes we are using both forward and reverse reads
    # So this is for reads that don't have mispriming issues
    
    fwdKey, revKey = sorted(readPairDict.keys())
    fwdSeq = readPairDict[fwdKey][0]
    revSeq = readPairDict[revKey][0]
    
    # Determine the 5' and 3' barcodes relative to each read
    # {read dir: {"end of sequence":[]}}
    allIDs = {"F":{"5":[], "3":[]}, "R":{"5":[], "3":[]}}
    allFoundBarcodes = {"F":{"5":[], "3":[]}, "R":{"5":[], "3":[]}} ###
    
    # Look for both barcodes in eithe read
    for direction, seq, allBarcodes in zip(["F", "R"], [fwdSeq, revSeq], \
                                [orderedFwdBarcodes, orderedRevBarcodes]):
        barcode5Region = seq[pred5BarcodeRange[0]:pred5BarcodeRange[1]]
        barcode3Region = seq[pred3BarcodeRange[0]:pred3BarcodeRange[1]]
        for i, degenerateBarcodeSet in enumerate(allBarcodes):
            found5 = False
            found3 = False
            for barcode in degenerateBarcodeSet:
                if not found5 and barcode in barcode5Region: 
                    allIDs[direction]["5"].append(i+1)
                    allFoundBarcodes[direction]["5"].append(barcode) ###
                    found5 = True
                if not found3 and barcode in barcode3Region:
                    allIDs[direction]["3"].append(i+1)
                    allFoundBarcodes[direction]["3"].append(barcode) ###
                    found3 = True
                if found5 and found3:
                    break # Both are matched, skip to the next barcode ID
            
    # Actually call what the barcodes were
    finalBarcodeCall = []
    all5IDs = allIDs["F"]["5"] + allIDs["R"]["3"]
    all3IDs = allIDs["F"]["3"] + allIDs["R"]["5"]
    
    for candidates in [all5IDs, all3IDs]:
        if len(candidates) != 0:
            possibleBarcode = most_common(candidates)
            if len(possibleBarcode) > 1: # more than 1 barcode is the most common
                finalBarcodeCall.append(None)
            else: # only one
                finalBarcodeCall.append(possibleBarcode[0])
        else: # No barcodes identified
            finalBarcodeCall.append(None)

    return(finalBarcodeCall) 

def barcode_identifier_JA20290_old(ReadPairDict, OrderedFwdBarcodes, OrderedRevBarcodes, \
                       pred_second_range_fwd = (190, 231), pred_second_range_rev = (195, 260), barcode_distance=238):
    # Identifies the forward and reverse barcodes by searching through the whole sequence
    # V2: takes best estimate of where the barcode is located in the Reverse and Forward reads
    # pred_second_range_fwd: min and max of range, where we look for the second barcode of fwd read
    # same for pred_second_range_rev, but for reverse read
    # Index is relative to 5' -> 3' directional indexing of reads
    # For JA20290, the reverse reads only catch 1 barcode
    
    # Set the bounds for where to look for 2nd barcode
    barcode_length = len(OrderedFwdBarcodes[0][0])
    
    fwd_key, rev_key = sorted(ReadPairDict.keys())
    fwd_seq = ReadPairDict[fwd_key][0]
    rev_seq = ReadPairDict[rev_key][0]
    
    fwdID = {"F":[], "R":[]} # which read showed the barcode 
    fwdID_idx = {"F":[], "R":[]} # position counting from the 5' end of the read that has barcode
    revID = {"F":[], "R":[]} 
    revID_idx = {"F":[], "R":[]} 
    
    # Forward and Reverse are treated differently now
    for direction, seq, all_barcodes in zip(["F", "R"],[fwd_seq, rev_seq], \
                                            [OrderedFwdBarcodes, OrderedRevBarcodes]):
        if direction == "F":
            seq_1 = seq[0:2*barcode_length] # This is now narrowed down to a specific region
        else:
            seq_1 = seq[pred_second_range_rev[0]: pred_second_range_rev[1]] 
        for i, barcode_set in enumerate(all_barcodes):
            for barcode in barcode_set: # All degenerate barcodes for one OG barcode
                if barcode in seq_1:
                    # locate what index it was found
                    if direction == "F": fwdID_idx[direction].append(seq_1.index(barcode))
                    else: fwdID_idx[direction].append(seq_1.index(barcode)+pred_second_range_rev[0])
                    
                    fwdID[direction].append(i+1)
                    break
                    
        # Only for the forward reads, attempt to find a second barcode
        if direction == "F":
            # Determine a range to look for possible 2nd barcodes
            if len(fwdID[direction]) != 0: # Have found fwd barcodes
                seq_2_range = [1000, 0] # Min max range to look at
                for idx in fwdID_idx[direction]:
                    lower_idx, upper_idx = idx, idx+barcode_length
                    second_barcode_lower_idx = barcode_distance + lower_idx + math.floor(barcode_length/2)
                    seq_2_range[0] = min(seq_2_range[0], second_barcode_lower_idx)
                    seq_2_range[1] = max(seq_2_range[1], second_barcode_lower_idx + 2*barcode_length)
            else: # didn't find first barcode
                seq_2_range = list(pred_second_range_fwd)
                
            seq_2 = seq[seq_2_range[0]:seq_2_range[1]]
            for i, barcode_set in enumerate(all_barcodes):
                for j, barcode in enumerate(barcode_set):
                    if barcode in seq_2: 
                        revID_idx[direction].append(seq_2.index(barcode) + seq_2_range[0])
                        revID[direction].append(i+1)
                        break
                        
                
    # Barcode-calling - determine what we're going to call the fwd and rev barcodes
    # If one is empty, default to the other
    # See if a common barcode exists; also the most commonly cited barcode
    # If there's a disagreement or neither determined it, call it an unknown barcode
    
    final_barcode_call = []
    for barcode_ID, barcode_idx in zip([fwdID, revID], [fwdID_idx, revID_idx]):
        candidates = barcode_ID["F"] + barcode_ID["R"]
        candidates = [x for x in candidates]
        if len(candidates) != 0:
            possible_bc = most_common(candidates)
            if len(possible_bc) > 1: # conflicting barcode calls; can't report
                final_barcode_call.append(None)
            else: # choices narrowed to one
                final_barcode_call.append(possible_bc[0])
        else: # all went to 0
            final_barcode_call.append(None)
    
    return(final_barcode_call)

def cell_bit_caller_old(fwd_seq, rev_seq, fwd_cells, rev_cells, fwd_bits_gapped, rev_bits_gapped):
    # Finds the value of the bits for each cell using results from pairwise2 alignment
    # fwd and rev cells are the unmutated constant regions for each cell
    # fwd and rev bits_gapped are the gapped variations on the 11-nt region with the bit nucleotide
    
    # Align fwd and rev using larger regions
    bit_id = []
    for i, (fwd_cell, rev_cell) in enumerate(zip(fwd_cells, rev_cells)):
        fwd_align = pairwise2.align.localms(fwd_cell, fwd_seq, 1, -0.5, -0.5, -0.5)
        _, _, fwd_score, _, _ = fwd_align[0]
        fwd_result_str = str(pairwise2.format_alignment(*fwd_align[0]))
        fwd_aligned_ref = fwd_result_str.split("\n")[0].split()[1]
        fwd_aligned_seq = fwd_result_str.split("\n")[2].split()[1]

        rev_align = pairwise2.align.localms(rev_cell, rev_seq, 1, -0.5, -0.5, -0.5)
        _, _, rev_score, _, _ = rev_align[0]
        rev_result_str = str(pairwise2.format_alignment(*rev_align[0]))
        rev_aligned_ref = rev_result_str.split("\n")[0].split()[1]
        rev_aligned_seq = rev_result_str.split("\n")[2].split()[1]
        
        # Compare the two scores, use it to narrow down where to look for bit values
        if fwd_score >= rev_score: # If they're the same, just pick fwd
            ref_seq = fwd_aligned_ref
            aligned_seq = fwd_aligned_seq
            selected_bits = fwd_bits_gapped[i]

        else: # Use the reverse
            ref_seq = rev_aligned_ref
            aligned_seq = rev_aligned_seq
            selected_bits = rev_bits_gapped[i]
        
        # Use list of gapped bit seqs to find bit values
        for bit in selected_bits:
            # Look for any instance of this bit region
            if bit in ref_seq: 
                start = ref_seq.index(bit)
                sample_bit_seq = aligned_seq[start: start+len(bit)]
                sample_bit_seq = "".join(sample_bit_seq.split("-"))
                sample_bit_nt = sample_bit_seq[math.floor(len(sample_bit_seq)/2)]
                ref_bit_seq = ref_seq[start: start+len(bit)]
                ref_bit_seq = "".join(ref_bit_seq.split("-"))
                ref_bit_nt = ref_bit_seq[math.floor(len(ref_bit_seq)/2)]
                
                # Check if they're the same
                if sample_bit_nt == ref_bit_nt:
                    bit_id.append(0)
                else: # Bits don't agree
                    # All "1" values ought to change to G for Fwd and C for Rev
                    if (fwd_score >= rev_score and sample_bit_nt == "G") or \
                        (fwd_score < rev_score and sample_bit_nt == "C"): # Rev
                        bit_id.append(1)
                    else: # Don't have the correct mutation; it's invalid
                        bit_id.append("?")
                
                break # Already found bit, shouldn't occur twice
        
        if len(bit_id) == i: # Didn't find a match to bit, this region probably didn't qualify the read
            print("Didn't find bit region for cell {}".format(i+1))
            bit_id.append("?")
        
    return(bit_id)

def count_correct_values_bitbybit(readsDictBySIMD, expectedBitsDict):
    # input: PlottedReads[nNGS]
    # expectedBitsDict = {1:"0001", 2:"0010",...} # keys are SIMD barcode indices, not values
    
    # Pad expectedBits if it doesn't expect anything for some barcodes
    if len(expectedBitsDict.keys()) < 16:
        for nSIMD in np.arange(1,16+1,1):
            try:
                expectedBitsDict[nSIMD];
            except KeyError:
                expectedBitsDict[nSIMD] = None
                
    possibleValues = ["0", "1", "?"]
    totalDigitCounts = {}
    for nSIMD in sorted(readsDictBySIMD.keys()):
        listOfReadBits = readsDictBySIMD[nSIMD]
        
        if len(listOfReadBits) == 0 or expectedBitsDict[nSIMD] == None: 
            # no reads for this init val
            totalDigitCounts[nSIMD] = None
            
        else:
            # Count up bits for each cell and normalize
            nReads = len(listOfReadBits)
            relativeDigitCount = [{x: 0 for x in possibleValues} for _ in range(4)]
            for readname, bits in listOfReadBits:
                for digitindex, value in enumerate(bits):
                    relativeDigitCount[digitindex][value] += 1
            relativeDigitCountNormalized = [{key: value/nReads for key, value in digit.items()} \
                                            for digit in relativeDigitCount]
            # Call the 4-bit value
            bitResult = ""
            for digit in relativeDigitCountNormalized:
                bitResult += max(digit, key=digit.get)
            totalDigitCounts[nSIMD] = bitResult
        
    # Determine if the correct values is present
    expectedBits = [expectedBitsDict[nSIMD] for nSIMD in sorted(readsDictBySIMD.keys())]
    actualBits = [totalDigitCounts[nSIMD] for nSIMD in sorted(readsDictBySIMD.keys())]
    matches = {}
    for x,y,nSIMD in zip(actualBits,expectedBits, sorted(readsDictBySIMD.keys())):
        if x == None or y == None: # no entry for this one
            matches[nSIMD] = None
        else:
            matches[nSIMD] = x == y
    
    return(matches, totalDigitCounts)

def count_correct_values_byread(readsDictBySIMD, expectedBitsDict):
    # input: PlottedReads[nNGS]
    # expectedBits: the correct 4-bit value you expect
    
    # Pad expectedBits if it doesn't expect anything for some barcodes
    if len(expectedBitsDict.keys()) < 16:
        for nSIMD in np.arange(1,16+1,1):
            try:
                expectedBitsDict[nSIMD];
            except KeyError:
                expectedBitsDict[nSIMD] = None
    
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    
    totalDigitCounts = {}
    for nSIMD in sorted(readsDictBySIMD.keys()):
        counts = [0 for _ in range(len(four_bits))]
        partial_bits = [] # Reads that had ? in them
        list_of_read_bits = readsDictBySIMD[nSIMD]
        
        if len(list_of_read_bits) == 0 or expectedBitsDict[nSIMD] == None: 
            # no reads for this init val
            totalDigitCounts[nSIMD] = None
            
        else:
            for readname, bit in list_of_read_bits:
                try:
                    index = four_bits.index(bit)
                    counts[index] += 1
                except ValueError: # not fully determined
                    partial_bits.append((readname, bit))
            # Call the most common read
            totalDigitCounts[nSIMD] = four_bits[counts.index(max(counts))]
    
    # Determine if the correct values is present
    expectedBits = [expectedBitsDict[nSIMD] for nSIMD in sorted(readsDictBySIMD.keys())]
    actualBits = [totalDigitCounts[nSIMD] for nSIMD in sorted(readsDictBySIMD.keys())]
    matches = {}
    for x,y,nSIMD in zip(actualBits,expectedBits, sorted(readsDictBySIMD.keys())):
        if x == None or y == None: # no entry for this one
            matches[nSIMD] = None
        else:
            matches[nSIMD] = x == y
  
    return(matches, totalDigitCounts)


def bitdistribution_matrix(readsDictBySIMD):
    # Makes a n x 16 matrix with the bit distribution (between 0 and 1) in it
    
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    orderednSIMDs = sorted(readsDictBySIMD.keys())
    
    # Bit distribution for control
    distroMatrix = []
    for nSIMD in orderednSIMDs:
        distroVector = [0 for _ in four_bits]
        for readname, bit in readsDictBySIMD[nSIMD]:
            try:
                index = four_bits.index(bit)
                distroVector[index] += 1
            except ValueError: # not fully determined
                continue 
        
        includednReads = sum(distroVector) # exclude partially determined
        try:
            distroMatrix.append([x/includednReads for x in distroVector])
        except ZeroDivisionError:
            distroMatrix.append([0 for _ in four_bits])
    
    return(np.array(distroMatrix))

def avg_matrix_TV(mat1,mat2):
    # Calculates the average total variation distance (1/2 of the L1 norm)
    # of two matrices
    # Has the interpretation of "Given any single read, what's the probability
    # it will produce the wrong answer?"
    # As long as this probably > 0.5, then given some reasonable number of NGS 
    # reads (100+), we can determine the correct answer ("amplification")
    
    dimensions = np.shape(mat1)
    TVs = []
    for i in range(dimensions[0]):
        TV = 0
        for j in range(dimensions[1]):
             TV += abs(mat1[i][j]-mat2[i][j])/2
        TVs.append(TV)
    return({"min":np.min(TVs), "mean":np.mean(TVs), "max":np.max(TVs)})

def make_report(plottedReads, indices, jobID, labels, idealAlgoMatrix, idealCtrlMatrix, title_y=1):
        
    for (i,j), label in zip(indices, labels):
        print("{} {}".format(jobID, label))

        rawCtrlMatrix = bitdistribution_matrix(plottedReads[i])
        rawCompMatrix = bitdistribution_matrix(plottedReads[j])
        orderednSIMDs = sorted(plottedReads[i].keys()) # assuming is the same for ctrl & comp
        
        # Pad matrices if not square
        padded = False
        if len(rawCtrlMatrix) != len(rawCtrlMatrix[0]):
            padded = True
            newCtrlMatrix = []
            newCompMatrix = []
            nrow = 0
            for r in range(16):
                if r+1 not in orderednSIMDs:
                    newCtrlMatrix.append(np.zeros(16))
                    newCompMatrix.append(np.zeros(16))
                else:
                    newCtrlMatrix.append(rawCtrlMatrix[nrow])
                    newCompMatrix.append(rawCompMatrix[nrow])
                    nrow += 1
            ctrlMatrix, compMatrix = newCtrlMatrix, newCompMatrix
        else:
            ctrlMatrix, compMatrix = rawCtrlMatrix, rawCompMatrix
        
        if not padded: # Will not be invertible otherwise; calculate transformation
            ctrlMatrix_inv = np.linalg.inv(ctrlMatrix)
            transMatrix = ctrlMatrix_inv @ compMatrix # Defined as data x transformation (transmat on right)
            simCompMatrix = ctrlMatrix @ idealAlgoMatrix # Simulated ideal comp given observed ctrl 
        
        if padded: # Reduce back down for report and visualization; do the same for ideal matrices
            finalCtrlMatrix = []
            finalCompMatrix = []
            finalIdealAlgoMatrix = []
            finalIdealCtrlMatrix = []
            for nSIMD in orderednSIMDs:
                finalCtrlMatrix.append(ctrlMatrix[nSIMD-1])
                finalCompMatrix.append(compMatrix[nSIMD-1])
                finalIdealAlgoMatrix.append(idealAlgoMatrix[nSIMD-1])
                finalIdealCtrlMatrix.append(idealCtrlMatrix[nSIMD-1])
        else:
            finalCtrlMatrix = ctrlMatrix
            finalCompMatrix = compMatrix
            finalIdealAlgoMatrix = idealAlgoMatrix
            finalIdealCtrlMatrix = idealCtrlMatrix
            
        # Calculate Total variation distance
        results = avg_matrix_TV(finalCompMatrix,finalIdealAlgoMatrix) # idealCount is also idealCtrl x idealCount 
        print("Mean Overall Error: {} (min={}, max={})".format(\
                    round(results['mean'],4), round(results['min'],4), round(results['max'],4)))

        results = avg_matrix_TV(finalCtrlMatrix,finalIdealCtrlMatrix)
        print("Mean Error in Control: {} (min={}, max={})".format(\
                    round(results['mean'],4), round(results['min'],4), round(results['max'],4)))
        
        if not padded:
            results = avg_matrix_TV(finalCompMatrix, simCompMatrix)
            print("Mean Error in Computation: {} (min={}, max={})".format(\
                        round(results['mean'],4), round(results['min'],4), round(results['max'],4)))

        # Make heatmaps
        matrix_heatmap(finalCompMatrix, finalIdealAlgoMatrix, title_adjust=title_y, \
                       title="Experimental Computation Matrix\n"+label, select_rows=[x-1 for x in orderednSIMDs])
        if not padded:
            matrix_heatmap(transMatrix, finalIdealAlgoMatrix, title_adjust=title_y, \
                           title="Simulated Computation from Ideal Control\n"+label, select_rows=[x-1 for x in orderednSIMDs])
            
        matrix_heatmap(finalCtrlMatrix, finalIdealCtrlMatrix, title_adjust=title_y, \
                       title="Experimental Control Matrix\n"+label, select_rows=[x-1 for x in orderednSIMDs])
        matrix_heatmap(finalIdealCtrlMatrix, finalIdealCtrlMatrix, title_adjust=title_y, \
                       title="Ideal Control Matrix\n"+label, select_rows=[x-1 for x in orderednSIMDs])
        matrix_heatmap(finalIdealAlgoMatrix, finalIdealAlgoMatrix, title_adjust=title_y, \
                       title="Ideal Computation Matrix\n"+label, select_rows=[x-1 for x in orderednSIMDs])
        if not padded:
            matrix_heatmap(simCompMatrix, finalIdealAlgoMatrix, title_adjust=title_y, \
                           title="Simulated Ideal Computation from Experimental Control\n"+label, select_rows=[x-1 for x in orderednSIMDs])

    return
    