from itertools import combinations, product
import numpy as np
import math
from Bio import pairwise2
from matplotlib.patches import Patch, Rectangle
import matplotlib.pyplot as plt
import matplotlib
import gzip
from matplotlib.gridspec import GridSpec 
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

def simd_ngs_heatmap_byexp(barcode_to_determined_bits, sample_index_to_barcode, \
                            initial_bits, expected_bits):
    # Like the Sanger result heatmaps, without digitization, arranged by experiment for easier viewing
    # Inputs are barcode mapping to determined bits, index mapping to barcode
    # Different from the previous as it doesn't necessarily consider both reads in a pair
    # The # of ? bits are indicated by size of rectangle
    # Updated 5/17/21: Can accomodate SIMD barcodes now, see "readname_bit" var
    
    possible_values = ["0", "1", "?"]
    n_samples = len(initial_bits)
    
    # Set up variables for later plotting
    actual_digit_count_normed = []
    unknown_count_normed = []
    reads_counted = []
    
    for i, barcode_pair in enumerate(sample_index_to_barcode):
        fwd_id, rev_id = barcode_pair

        # Organize data into list of dicts
        n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
        reads_counted.append(n_reads)
        if n_reads == 0:
            barcode_digit_count_normed = [np.nan for _ in range(4)]
            barcode_unknown_count_normed = [0 for _ in range(4)]
        else:
            rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
            for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                readname, bit = readname_bit[0], readname_bit[1]
                for digit_index, value in enumerate(bit):
                    rel_digit_count[digit_index][value] += 1
            norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
            rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]
            
            barcode_digit_count_normed = []
            for v in rel_digit_count_normed:
                try:
                    barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                except ZeroDivisionError:
                    barcode_digit_count_normed.append(np.nan)
            barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]
        
        actual_digit_count_normed.append(barcode_digit_count_normed)
        unknown_count_normed.append(barcode_unknown_count_normed)
        
    # Prepare initial and expected values
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [initial_values, expected_values, actual_digit_count_normed]
    all_labels = ["Initial", "Expected", "Actual"]
    
    # Plot
    plt.close()
    fig, axes = plt.subplots(1, 3, figsize=(15, n_samples))
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    for i in range(3):
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = axes[i]

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
    # Turn on label for just first plot
    axes[0].set_yticklabels(initial_bits, minor=True)
       
    # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction  
    for r, row in enumerate(unknown_count_normed):
        for c, value in enumerate(row):
            if value > 0:
                half_edge = math.sqrt(value)/2
                rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                axes[-1].add_patch(rect)
    
    # Add right hand labels to Actual as read numbers
    axes[-1].set_yticklabels(["n = {}".format(x) for x in reads_counted], minor=True)
    axes[-1].yaxis.tick_right()
    
    # Make heatmap legend
    fig.subplots_adjust(right=0.76)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
    cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
    cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.show()
    
    return

def simd_ngs_heatmap_initial_actual(barcode_to_determined_bits, sample_index_to_barcode, \
                            expected_bits, legend=True, plot_adjust=0.76, title=None, \
                            title_y=1.02):
    # Simplest version of the viridis heatmap, just takes expected and 
    # the experimental outcome
    # Note that these plots take only both-NGS barcode reads
    # For presentation purposes. Can toggle legend on/off
    # The # of ? bits are indicated by size of rectangle
    
    possible_values = ["0", "1", "?"]
    n_samples = len(expected_bits)
    
    # Set up variables for later plotting
    actual_digit_count_normed = []
    unknown_count_normed = []
    reads_counted = []
    
    for i, barcode_pair in enumerate(sample_index_to_barcode):
        fwd_id, rev_id = barcode_pair

        # Organize data into list of dicts
        n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
        reads_counted.append(n_reads)
        if n_reads == 0:
            barcode_digit_count_normed = [np.nan for _ in range(4)]
            barcode_unknown_count_normed = [0 for _ in range(4)]
        else:
            rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
            for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                readname, bit = readname_bit[0], readname_bit[1]
                for digit_index, value in enumerate(bit):
                    rel_digit_count[digit_index][value] += 1
            norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
            rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]
            
            barcode_digit_count_normed = []
            for v in rel_digit_count_normed:
                try:
                    barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                except ZeroDivisionError:
                    barcode_digit_count_normed.append(np.nan)
            barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]
        
        actual_digit_count_normed.append(barcode_digit_count_normed)
        unknown_count_normed.append(barcode_unknown_count_normed)
        
    # Prepare initial and expected values
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [expected_values, actual_digit_count_normed]
    all_labels = ["Expected", "Actual"]
    
    # Plot
    plt.close()
    fig, axes = plt.subplots(1, 2, figsize=(9, n_samples))
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    for i in range(2):
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = axes[i]

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
    # Turn on label for just first plot
    axes[0].set_yticklabels(expected_bits, minor=True)
       
    # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction  
    for r, row in enumerate(unknown_count_normed):
        for c, value in enumerate(row):
            if value > 0:
                half_edge = math.sqrt(value)/2
                rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                axes[-1].add_patch(rect)
    
    # Add right hand labels to Actual as read numbers
    axes[-1].set_yticklabels(["n = {}".format(x) for x in reads_counted], minor=True)
    axes[-1].yaxis.tick_right()
    
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(right=plot_adjust)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
        cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    if title != None:
        plt.suptitle(title, y=title_y)
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.show()
    
    return

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

def sample_to_bits_breakdown(list_of_read_bits, label, expected_bits=None):
    # Given a list of (readname, determined bit value), plots distribution of all bit values
    # For now, leaves out ? reads
    # This is different from the other sample_to_bits_breakdown functions in that
    # it only expects to show results for a single sample

    # Tally up all the bits
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    counts = [0 for _ in range(len(four_bits))]
    partial_bits = [] # Reads that had ? in them
    for readname, bit in list_of_read_bits:
        try:
            index = four_bits.index(bit)
            counts[index] += 1
        except ValueError: # not fully determined
            partial_bits.append((readname, bit))
            
    normedcounts = [100*x/len(list_of_read_bits) for x in counts]

    # Plot as bar plot
    plt.close()
    fig, ax = plt.subplots(1,1, figsize=(len(four_bits),6))
    width = 0.75
    x = np.arange(0, len(four_bits))
    ax.bar(x, normedcounts, width, color="#00a9b7")
    for i, v in enumerate(normedcounts):
        ax.text(i, v+max(normedcounts)*0.025, str(round(v,2)), color='black', fontweight='bold', ha="center")
    ax.set_xticks(x)
    ax.set_xticklabels(four_bits)
    ax.set_title(label+" (total n={}, fully decoded={})".format(len(list_of_read_bits), sum(counts)), y=1.02)
    ax.set_ylim(0, 110)
    ax.set_yticks(np.arange(0,110,20))
    ax.set_yticklabels(np.arange(0,110,20))
    ax.set_ylabel("Percent of all reads (per barcoded sample)")
    if expected_bits != None:
        if len(expected_bits) == 2: # We will assume the first is the control, second is the +1
            ax.get_xticklabels()[expected_bits[0]-1].set_color("#f8971f") # ctrl is orange
            ax.get_xticklabels()[expected_bits[0]-1].set_weight("bold") 
            ax.get_xticklabels()[expected_bits[1]-1].set_color("#a6cd57") # +1 is green
            ax.get_xticklabels()[expected_bits[1]-1].set_weight("bold") 
            legend_handles = [Patch(color="#f8971f"), Patch(color="#a6cd57")]
            ax.legend(legend_handles, ["Ctrl", "+1"], fontsize=16, loc=1, ncol=2, framealpha=0.5)
        elif len(expected_bits) == 3: # We'll assume first is ctrl, second +1, and third +1+1
            ax.get_xticklabels()[expected_bits[0]-1].set_color("#f8971f") # ctrl is orange
            ax.get_xticklabels()[expected_bits[0]-1].set_weight("bold") 
            ax.get_xticklabels()[expected_bits[1]-1].set_color("#a6cd57") # +1 is green
            ax.get_xticklabels()[expected_bits[1]-1].set_weight("bold") 
            ax.get_xticklabels()[expected_bits[2]-1].set_color("#ffd600") # +1+1 is yellow
            ax.get_xticklabels()[expected_bits[2]-1].set_weight("bold") 
            legend_handles = [Patch(color="#f8971f"), Patch(color="#a6cd57"), Patch(color="#ffd600")]
            ax.legend(legend_handles, ["Ctrl", "+1", "+1+1"], fontsize=16, loc=1, ncol=3, framealpha=0.5)
        else: # No assumption about what is what, so make everything red-bold
            for j in expected_bits:
                ax.get_xticklabels()[j-1].set_color("red")
                ax.get_xticklabels()[j-1].set_weight("bold")
    ax.grid(True)
    plt.show()
    return(counts, len(partial_bits))
    
def sample_to_bits_breakdown_ctrl_expt(list_of_read_bits_ctrl, list_of_read_bits_expt, \
    label, expected_ctrl=False, expected_expt=False, bar_labels=["Ctrl", "Expt"]):
    # Given a list of (readname, determined bit value), plots distribution of all bit values
    # For now, leaves out ? reads
    # This is different from the other sample_to_bits_breakdown functions in that
    # it only expects to show results for two samples (ctrl and corresponding expt)
    # UPDATE 7/18/21: More generalized title and labels (bar_labels)

    # Tally up all the bits
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    counts_ctrl = [0 for _ in range(len(four_bits))]
    partial_bits_ctrl = []
    for readname, bit in list_of_read_bits_ctrl:
        try:
            index = four_bits.index(bit)
            counts_ctrl[index] += 1
        except ValueError: # not fully determined
            partial_bits_ctrl.append((readname, bit))
    norm_ctrl = sum(counts_ctrl)
    if norm_ctrl > 0:
        counts_ctrl = [100*x/norm_ctrl for x in counts_ctrl]
    else: # No data
        counts_ctrl = [0 for x in counts_ctrl]
            
    counts_expt = [0 for _ in range(len(four_bits))]
    partial_bits_expt = []
    for readname, bit in list_of_read_bits_expt:
        try:
            index = four_bits.index(bit)
            counts_expt[index] += 1
        except ValueError: # not fully determined
            partial_bits_expt.append((readname, bit))
    norm_expt = sum(counts_expt)
    counts_expt = [100*x/norm_expt for x in counts_expt]
    
    # Plot as bar plot
    plt.close()
    colors = ["#f8971f", "#ffd600", "#a6cd57", "#579d42", "#00a9b7"]
    fig, ax = plt.subplots(1,1, figsize=(len(four_bits),6))
    width = 0.7
    gap = 0.1
    x = np.arange(0, len(four_bits)*2, 2)
    ax.bar(x, counts_ctrl, width, color=colors[0])
    ax.bar(x+width+gap, counts_expt, width, color=colors[2])
    for i, v in enumerate(counts_ctrl):
        ax.text(i*2, v+max(counts_ctrl)*0.025, str(round(v,2)), color='black', fontweight='bold', ha="center", fontsize=12, rotation=45)
        w = counts_expt[i]
        ax.text(i*2+width+gap, w+max(counts_expt)*0.025, str(round(w, 2)), color='red', fontweight='bold', ha="center", fontsize=12, rotation=45)
    ax.set_xticks(x)
    ax.set_xticklabels(four_bits)
    ax.set_ylabel("Percent of all reads (per barcoded sample)")
    ax.set_title(label+"\n({} n={}, {} n={})".format(bar_labels[0],len(list_of_read_bits_ctrl), \
                                                    bar_labels[1],len(list_of_read_bits_expt)), y=1.02)
    
    legend_handles = [Patch(color=colors[0]), Patch(color=colors[2])]
    ax.legend(legend_handles, bar_labels, fontsize=16, loc=1, ncol=n_bits+1, framealpha=0.5)
    #ax.set_ylim(0, max(max(counts_ctrl), max(counts_expt))*1.1)
    ax.set_ylim(0,110)
    ax.set_yticks(np.arange(0,110,20))
    ax.set_yticklabels(np.arange(0,110,20))
    
    if expected_ctrl != False:
        ax.get_xticklabels()[four_bits.index(expected_ctrl)].set_color("#f8971f")
        ax.get_xticklabels()[four_bits.index(expected_ctrl)].set_weight("bold")
    
    if expected_expt != False:
        ax.get_xticklabels()[four_bits.index(expected_expt)].set_color("#a6cd57")
        ax.get_xticklabels()[four_bits.index(expected_expt)].set_weight("bold")

    ax.grid(True)
    plt.show()
         
    return(counts_ctrl, counts_expt)

def sample_to_bits_breakdown_ctrl_expt_2(list_of_read_bits_ctrl, \
                                        list_of_read_bits_expt, \
                                        list_of_read_bits_expt_2, \
                                        label, expected_ctrl=False,
                                        expected_expt=False, expected_expt_2=False,\
                                        bar_labels=["Ctrl", "Expt +1", "Expt +1+1"]):
    # Given a list of (readname, determined bit value), plots distribution of all bit values
    # For now, leaves out ? reads
    # This version includes 2 sets of expt data - for comparing +1 and +1+1
    # Update 7/18/21: Now has a more generalized title text
    
    # Tally up all the bits
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    counts_ctrl = [0 for _ in range(len(four_bits))]
    
    # Control
    partial_bits_ctrl = []
    for readname, bit in list_of_read_bits_ctrl:
        try:
            index = four_bits.index(bit)
            counts_ctrl[index] += 1
        except ValueError: # not fully determined
            partial_bits_ctrl.append((readname, bit))
    norm_ctrl = sum(counts_ctrl)
    counts_ctrl = [100*x/norm_ctrl for x in counts_ctrl]
    
    # Experiment (+1)
    counts_expt = [0 for _ in range(len(four_bits))]
    partial_bits_expt = []
    for readname, bit in list_of_read_bits_expt:
        try:
            index = four_bits.index(bit)
            counts_expt[index] += 1
        except ValueError: # not fully determined
            partial_bits_expt.append((readname, bit))
    norm_expt = sum(counts_expt)
    counts_expt = [100*x/norm_expt for x in counts_expt]
    
    # Experiment 2 (+1 +1)
    counts_expt_2 = [0 for _ in range(len(four_bits))]
    partial_bits_expt_2 = []
    for readname, bit in list_of_read_bits_expt_2:
        try:
            index = four_bits.index(bit)
            counts_expt_2[index] += 1
        except ValueError: # not fully determined
            partial_bits_expt_2.append((readname, bit))
    norm_expt_2 = sum(counts_expt_2)
    counts_expt_2 = [100*x/norm_expt_2 for x in counts_expt_2]
    
    # Plot as bar plot
    plt.close()
    colors = ["#a6cd57", "#00a9b7","#f8971f","#579d42","#005f86","#ffd600"]
    fig, ax = plt.subplots(1,1, figsize=(len(four_bits),6))
    width = 0.7
    gap = 0.1
    x = np.arange(0, len(four_bits)*3, 3)
    ax.bar(x, counts_ctrl, width, color=colors[0])
    ax.bar(x+width+gap, counts_expt, width, color=colors[2])
    ax.bar(x+(width+gap)*2, counts_expt_2, width, color=colors[1])
    
    # Label text on top of bar plots
    text_colors = ["#e39220", "#8db34f", "#d6ba03"]
    all_counts = [counts_ctrl, counts_expt, counts_expt_2]
    for i, xpos in enumerate(x):
        for j, (color, counts) in enumerate(zip(text_colors, all_counts)):
            v = counts[i]
            ax.text(xpos+j*(width+gap), v+max(counts_ctrl)*0.025, str(round(v,2)), \
                    color=color, fontweight='bold', ha="center", fontsize=12, rotation=60)
        
    ax.set_xticks(x+width+gap)
    ax.set_xticklabels(four_bits)
    ax.set_ylabel("Percent of all reads (per barcoded sample)")
    ax.set_title(label+"\n({} n={}, {} n={}, {} n={})".format(\
                                                        bar_labels[0],\
                                                        len(list_of_read_bits_ctrl), \
                                                        bar_labels[1],\
                                                        len(list_of_read_bits_expt), \
                                                        bar_labels[2],
                                                        len(list_of_read_bits_expt_2)), y=1.02)
    
    legend_handles = [Patch(color=colors[0]), Patch(color=colors[2]), Patch(color=colors[1])]
    ax.legend(legend_handles, bar_labels, fontsize=16, loc=1, \
              ncol=n_bits+1, framealpha=0.5)
    #ax.set_ylim(0, max(max(counts_ctrl), max(counts_expt))*1.1)
    ax.set_ylim(0,110)
    ax.set_yticks(np.arange(0,110,20))
    ax.set_yticklabels(np.arange(0,110,20))
    
    if expected_ctrl != False:
        ax.get_xticklabels()[four_bits.index(expected_ctrl)].set_color("#e39220")
        ax.get_xticklabels()[four_bits.index(expected_ctrl)].set_weight("bold")
    
    if expected_expt != False:
        ax.get_xticklabels()[four_bits.index(expected_expt)].set_color("#8db34f")
        ax.get_xticklabels()[four_bits.index(expected_expt)].set_weight("bold")
        
    if expected_expt_2 != False:
        ax.get_xticklabels()[four_bits.index(expected_expt_2)].set_color("#d6ba03")
        ax.get_xticklabels()[four_bits.index(expected_expt_2)].set_weight("bold")

    ax.grid(True)
    plt.show()
         
    return(counts_ctrl, counts_expt, counts_expt_2)

def simd_ngs_heatmap_byexp_with_control(barcode_to_determined_bits, ctrl_sample_index_to_barcode, \
                            expt_sample_index_to_barcode, initial_bits, expected_bits):
    # Like the Sanger result heatmaps, without digitization, arranged by experiment for easier viewing
    # Inputs are barcode mapping to determined bits, index mapping to barcode
    # Different from the previous as it doesn't necessarily consider both reads in a pair
    # The # of ? bits are indicated by size of rectangle
    # UPDATE 10/19/20: Now includes control heatmap in addition to experiment heatmap
    # UPDATE 5/15/20: Now accounts for reads with SIMD barcode, which adds to the readname, bit tuple
    
    possible_values = ["0", "1", "?"]
    n_samples = len(initial_bits)
    
    # Set up variables for later plotting: Control bits
    ctrl_digit_count_normed = [] # was once actual_digit_count_normed
    ctrl_unknown_count_normed = []
    ctrl_reads_counted = []
    
    for i, barcode_pair in enumerate(ctrl_sample_index_to_barcode):
        fwd_id, rev_id = barcode_pair

        # Organize data into list of dicts
        n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
        ctrl_reads_counted.append(n_reads)
        if n_reads == 0:
            barcode_digit_count_normed = [np.nan for _ in range(4)]
            barcode_unknown_count_normed = [0 for _ in range(4)]
        else:
            rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
            for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                readname, bit = readname_bit[0], readname_bit[1] # Don't care about SIMD ID
                for digit_index, value in enumerate(bit):
                    rel_digit_count[digit_index][value] += 1
            norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
            rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]
            
            barcode_digit_count_normed = []
            for v in rel_digit_count_normed:
                try:
                    barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                except ZeroDivisionError:
                    barcode_digit_count_normed.append(np.nan)
            barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]
        
        ctrl_digit_count_normed.append(barcode_digit_count_normed)
        ctrl_unknown_count_normed.append(barcode_unknown_count_normed)
        
    # Set up variables for later plotting: Experiment bits 
    expt_digit_count_normed = [] # was once actual_digit_count_normed
    expt_unknown_count_normed = []
    expt_reads_counted = []
    
    for i, barcode_pair in enumerate(expt_sample_index_to_barcode):
        fwd_id, rev_id = barcode_pair

        # Organize data into list of dicts
        n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
        expt_reads_counted.append(n_reads)
        if n_reads == 0:
            barcode_digit_count_normed = [np.nan for _ in range(4)]
            barcode_unknown_count_normed = [0 for _ in range(4)]
        else:
            rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
            for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                readname, bit = readname_bit[0], readname_bit[1]
                for digit_index, value in enumerate(bit):
                    rel_digit_count[digit_index][value] += 1
            norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
            rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]
            
            barcode_digit_count_normed = []
            for v in rel_digit_count_normed:
                try:
                    barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                except ZeroDivisionError:
                    barcode_digit_count_normed.append(np.nan)
            barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]
        
        expt_digit_count_normed.append(barcode_digit_count_normed)
        expt_unknown_count_normed.append(barcode_unknown_count_normed)
        
    # Prepare initial and expected values
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [initial_values, ctrl_digit_count_normed, expected_values,  expt_digit_count_normed]
    all_labels = ["Initial", "Control", "Expected", "Experiment"]
    
    # Plot
    plt.close()
    fig = plt.figure(figsize=(15, n_samples), constrained_layout=True)
    outer_spec = fig.add_gridspec(1,2,wspace=0.1, width_ratios = [1,1])
    inner_spec_ctrl = outer_spec[0].subgridspec(1,2)
    inner_spec_expt = outer_spec[1].subgridspec(1,2)
    ax0 = fig.add_subplot(inner_spec_ctrl[0]) # initial values
    ax1 = fig.add_subplot(inner_spec_ctrl[1]) # ctrl
    ax2 = fig.add_subplot(inner_spec_expt[0]) # expected + 1
    ax3 = fig.add_subplot(inner_spec_expt[1]) # experiment
    axes = [ax0, ax1, ax2, ax3]#, ax4]
    
    #current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap = copy.copy(matplotlib.cm.get_cmap("viridis"))
    current_cmap.set_bad(color="gray")
    for i in range(4):
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = axes[i]

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
    # Turn on label for first plot and third plot
    axes[0].set_yticklabels(initial_bits, minor=True)
    axes[2].set_yticklabels(expected_bits, minor=True)
    
    # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction  
    for r, row in enumerate(ctrl_unknown_count_normed):
        for c, value in enumerate(row):
            if value > 0:
                half_edge = math.sqrt(value)/2
                rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                axes[1].add_patch(rect)
                
    for r, row in enumerate(expt_unknown_count_normed):
        for c, value in enumerate(row):
            if value > 0:
                half_edge = math.sqrt(value)/2
                rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                axes[3].add_patch(rect)
    
    # Add right hand labels to Control and Experiment as read numbers
    axes[1].set_yticklabels(["n = {}".format(x) for x in ctrl_reads_counted], minor=True, fontsize=14)
    axes[1].yaxis.tick_right() ### Add adjustment to accomodate different expt conditions
    
    axes[3].set_yticklabels(["n = {}".format(x) for x in expt_reads_counted], minor=True, fontsize=14)
    axes[3].yaxis.tick_right()
    
    # Make heatmap legend
    cax = fig.add_axes([axes[3].get_position().x1+0.11,axes[3].get_position().y0,0.02,axes[3].get_position().height])
    cbar_handle = fig.colorbar(handle, cax=cax, ticks=[0, 0.5, 1])
    cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
    cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
 
    plt.show()
    
    return

def sample_to_bits_breakdown_multi(list_of_lists_of_read_bits, \
                                    title_label, bar_labels, \
                                   expt_values=False, legend_loc=1):
    # Distribution of bits for each sample, for up to 6 different
    # samples
    # Everything is in the form of lists
    # Also reports partial bits in the title
    # expt_values is the manual option to highlight expected bit value
    # in bold and color; must be list
    
    # Tally up all the bits
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    nSamples = len(list_of_lists_of_read_bits)
    
    # Save all counts as a list of lists
    counts_all_percent = []
    partial_counts_all = []
    counts_all = []
    for list_of_read_bits in list_of_lists_of_read_bits:
        counts_this_sample = [0 for _ in range(len(four_bits))]
        partial_bits_this_sample = []
        for readname, bit in list_of_read_bits:
            try: 
                index = four_bits.index(bit)
                counts_this_sample[index] += 1
            except ValueError: # not fully determined; not included
                partial_bits_this_sample.append((readname, bit))
        norm_sample = sum(counts_this_sample)
        counts_all.append(norm_sample)
        counts_all_percent.append([100*x/norm_sample for x in counts_this_sample])
        partial_counts_all.append(len(partial_bits_this_sample))
    
    # Plot as bar plot
    plt.close()
    colors = ["#a6cd57", "#00a9b7","#f8971f","#579d42","#005f86","#ffd600"]
    fig, ax = plt.subplots(1,1, figsize=(len(four_bits),6))
    width = 0.7
    gap = 0.1
    ceiling = max([max(counts_all_percent[x]) for x in range(nSamples)])
    for i in range(nSamples):
        x = np.arange((width+gap)*i+0, 16*nSamples+(width+gap)*i, nSamples)
        ax.bar(x, counts_all_percent[i], color=colors[i])

        # Label text on top of bar plots
        for j, xpos in enumerate(x):
            value = counts_all_percent[i][j]
            ax.text(xpos, value+ceiling*0.025, str(round(value,2)),\
                   color=colors[i], fontweight="bold", ha="center", \
                    fontsize=12, rotation=60)
            
    ax.set_xticks(np.arange(0,16*nSamples,nSamples)+(width+gap)*(nSamples-1)/2)
    ax.set_xticklabels(four_bits)
    ax.set_ylabel("Percent of fully-decoded reads in sample")
    titleReport = ["{} n={} ({})".format(name, c, p) \
                   for name, c, p in \
                   zip(bar_labels, counts_all, partial_counts_all)]
    ax.set_title(title_label, y=1.02)

    legend_handles = [Patch(color=colors[x]) for x in range(nSamples)]
    ax.legend(legend_handles, titleReport, fontsize=16, loc=legend_loc, \
              ncol=n_bits+1, framealpha=0.5)
    ax.set_ylim(0,110)
    ax.set_yticks(np.arange(0,110,20))
    ax.set_yticklabels(np.arange(0,110,20))
    
    # Highlight "correct" values; only does up to 2, manually specify
    if expt_values != False:
        for i, (value, color) in enumerate(zip(expt_values, colors)):
            ax.get_xticklabels()[four_bits.index(value)].set_color(color)
            ax.get_xticklabels()[four_bits.index(value)].set_weight("bold")

    ax.grid(True)
    plt.show()
         
    return(counts_all, partial_counts_all)

def sample_to_bits_breakdown_multi_grid(list_of_lists_of_lists_of_read_bits, \
                                    list_of_title_label, bar_labels, suptitle, \
                                    dimension=(4,2), list_of_expt_values=False):
    # Same as sample_to_bits_breakdown_multi but plots several results together
    # Plots first down the column, then across the row
    
    # Tally up all the bits
    four_bits = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    
    plt.close()
    figRows, figCols = dimension
    fig, axes = plt.subplots(figRows, figCols, figsize=(figCols*12,figRows*5))
    colors = ["#a6cd57", "#00a9b7","#f8971f","#579d42","#005f86","#ffd600"]
    width = 0.7
    gap = 0.1
    
    for i, singlePlotData in enumerate(zip(list_of_lists_of_lists_of_read_bits, \
                                    list_of_title_label, list_of_expt_values)):
        list_of_lists_of_read_bits, title_label, expt_values = singlePlotData
        nSamples = len(list_of_lists_of_read_bits)
        if figRows > 1 and figCols > 1:
            ax = axes[i%figRows][math.floor(i/figRows)]
        else:
            ax = axes[i]
            
        # Save all counts as a list of lists
        counts_all_percent = []
        partial_counts_all = []
        counts_all = []
        for list_of_read_bits in list_of_lists_of_read_bits:
            counts_this_sample = [0 for _ in range(len(four_bits))]
            partial_bits_this_sample = []
            for readname, bit in list_of_read_bits:
                try: 
                    index = four_bits.index(bit)
                    counts_this_sample[index] += 1
                except ValueError: # not fully determined; not included
                    partial_bits_this_sample.append((readname, bit))
            norm_sample = sum(counts_this_sample)
            counts_all.append(norm_sample)
            counts_all_percent.append([100*x/norm_sample for x in counts_this_sample])
            partial_counts_all.append(len(partial_bits_this_sample))
    
        # Plot as bar plot
        ceiling = max([max(counts_all_percent[x]) for x in range(nSamples)])
        for i in range(nSamples):
            x = np.arange((width+gap)*i+0, 16*nSamples+(width+gap)*i, nSamples)
            ax.bar(x, counts_all_percent[i], color=colors[i]) 

            # Label text on top of bar plots
            for j, xpos in enumerate(x):
                value = counts_all_percent[i][j]
                ax.text(xpos, value+ceiling*0.025, str(round(value,2)),\
                       color=colors[i], fontweight="bold", ha="center", \
                        fontsize=12, rotation=60)

        ax.set_xticks(np.arange(0,16*nSamples,nSamples)+(width+gap)*(nSamples-1)/2)
        ax.set_xticklabels(four_bits, fontsize=12)
        ax.set_ylabel("Percent reads in sample")
        reportCounts = ["{} n={} ({})".format(name, c, p) \
                       for name, c, p in \
                       zip(bar_labels, counts_all, partial_counts_all)]
        ax.set_title(title_label, y=1.02) #+", Fully-decoded (partially-decoded)"

        legend_handles = [Patch(color=colors[x]) for x in range(nSamples)]
        ax.legend(legend_handles, reportCounts, fontsize=10, loc=1, \
                  ncol=n_bits+1, framealpha=0.5)
        ax.set_ylim(0,120)
        ax.set_yticks(np.arange(0,110,20))
        ax.set_yticklabels(np.arange(0,110,20))

        # Highlight "correct" values; only does up to 2, manually specify
        if expt_values != False:
            for i, (value, color) in enumerate(zip(expt_values,colors)):
                ax.get_xticklabels()[four_bits.index(value)].set_color(color)
                ax.get_xticklabels()[four_bits.index(value)].set_weight("bold")

        ax.grid(True)
    plt.suptitle(suptitle, fontsize=24, y=0.92+0.02*(figCols-1))
    plt.subplots_adjust(hspace=0.3, wspace=0.15)
    plt.show()
         
    return


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

def expected_error_by_value(barcode_to_determined_bits, data_indices, expected_bits, title,\
                            toss_undetermined=True, n_samples=8, n_bits=4):
    # Plots distribution of bit values (correlated 4-bit values) by error, reports expected error
    
    plt.close()
    plt.figure(figsize=(n_samples*2, 5))
    colors = ["#f8971f", "#ffd600", "#a6cd57", "#579d42", "#00a9b7"]
    labels = ["{} errors".format(k) for k in np.arange(n_bits+1)]
    label_x_pos = np.arange((n_bits+2)/2,n_samples*(n_bits+2), (n_bits+2))-1
    
    all_expected_error = []
    
    for i in range(n_samples):
        correct_value = expected_bits[i]
        fwd_idx, rev_idx = data_indices[i]
        bits_list = [x[1] for x in barcode_to_determined_bits[fwd_idx][rev_idx]]
        sorted_bits = {x:[] for x in range(n_bits+1)}
        for bit in bits_list:
            if toss_undetermined and "?" in bit:
                continue
            else:
                error = hamming_dist(bit, correct_value)
                sorted_bits[error] += bit
    
        data = [len(sorted_bits[x]) for x in range(n_bits+1)]
        total_n = sum(data)
        if total_n == 0: 
            normed_data = [None for _ in range(n_bits+1)] 
            print("Sample {} has 0 data!".format(i))
        else:
            normed_data = [x/total_n for x in data]
        
        # Plot data
        if total_n != 0: # Actual data
            width = 0.75
            for j in range(n_bits+1):
                x = i*(n_bits+2)+j
                plt.bar(x, normed_data[j], width, color=colors[j%5], label=labels[j])

            expected_error = np.sum([n*x for n,x in enumerate(normed_data)])
            plt.text(i*(n_bits+2), max(normed_data)*1.05, round(expected_error,2), \
                     color='black', ha="center")
            all_expected_error.append(expected_error)
        else: # Don't have actual data, so just plot a big black bar at 1.0
            width = n_bits/2.0
            x = i*(n_bits+2) + 2 
            plt.bar(x, 1, width, color="gray", align="center")
            plt.text(x, 1.05, "No data", color="black", ha="center")
            all_expected_error.append(None)
    
    legend_handles = []
    for j in range(n_bits+1):
        legend_handles.append(Patch(color=colors[j])) 
    legend_handles.append(Patch(color="gray"))
    plt.legend(legend_handles, labels+["No data"], fontsize=14, loc=1, ncol=n_bits+2)
    
    plt.xlim(-2, n_samples*(n_bits+2))
    plt.ylim(0,1.4)
    plt.title("{} Bit sequence error (Hamming distance) by expected value".format(title), y=1.05)
    plt.ylabel("Proportion of reads (Normalized)")
    plt.xlabel("Expected bit value")
    plt.xticks(label_x_pos, expected_bits, ha="center")
    plt.grid(True)
    plt.show()
        
    return(all_expected_error)

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

def simd_sanger_barplot(actual_controls, exp_names, legend=True):
    # Looks at one set of data (control or +1) only; presents it in bargraph format
    
    # Process actual data    
    actual_control_fwd, actual_control_rev = actual_controls
    
    # Controls
    def normalize_bits(results):
        results_normed = []
        for sample in results:
            sample_bits = []
            for x,y in sample:
                if x == "?" or y == "?":
                    sample_bits.append((x, y))
                else:
                    sample_bits.append((x/(x+y), y/(x+y)))
            results_normed.append(sample_bits)
        return(results_normed)
                
    actual_control_fwd_normed = normalize_bits(actual_control_fwd)
    actual_control_rev_normed = normalize_bits(actual_control_rev)
    
    # Too messy to do by list comp
    def combined_normalized_bits(fwd_normed, rev_normed):
        actual_normed = []
        for sample_fwd, sample_rev in zip(fwd_normed, rev_normed):
            sample_normed_fwd = []
            sample_normed_rev = []

            for fwd, rev in zip(sample_fwd, sample_rev):
                # Save values based on which (Fwd or Rev) is viable
                # Odd rows will be fwd, even will be rev
                if "?" not in fwd:
                    sample_normed_fwd.append((fwd[1]))
                else:
                    sample_normed_fwd.append(float('NaN'))
                    
                if "?" not in rev:
                    sample_normed_rev.append(rev[1])
                else:
                    sample_normed_rev.append(float('NaN'))

            actual_normed.append((sample_normed_fwd, sample_normed_rev))
        return(actual_normed)
    
    actual_control_normed = combined_normalized_bits(actual_control_fwd_normed, actual_control_rev_normed)
    
    # Plot
    width = 0.25
    n_samples = len(actual_controls)
    for i, (pt, name) in enumerate(zip(actual_control_normed, exp_names)):
        x_fwd, x_rev = pt
        
        plt.close()
        plt.figure(figsize=(8, 3))
        for pos, x in zip(np.arange(4), x_fwd):
            if not np.isnan(x):
                plt.bar(pos, x, width, color="Blue", label="Fwd Read")
            else:
                plt.bar(pos, 0.5, width, color="Gray", label="Invalid Fwd")
                
        for pos, x in zip(np.arange(4), x_rev):
            if not np.isnan(x):
                plt.bar(pos+1.25*width, x, width, color="Red", label="Red Read")
            else:
                plt.bar(pos+1.25*width, 0.5, width, color="Black", label="Invalid Rev")
        
        if legend:
            legend_handles = [Patch(color="Blue"), Patch(color="Gray"),\
                              Patch(color="Red"), Patch(color="Black")]
            legend_labels = ["Fwd Read", "Invalid Fwd", "Rev Read", "Invalid Rev"]
            plt.legend(legend_handles, legend_labels, bbox_to_anchor=(1,0.5), \
                       loc=6, edgecolor="black", title="Sanger")        
        plt.axhline(0.5, color="gray", linestyle="--")
        plt.grid(True)
        plt.ylim(0,1)
        plt.xlabel("Bit")
        plt.ylabel("Normalized bit value")
        plt.xticks(np.arange(4)+0.625*width, np.arange(1,5))
        plt.title(name)
        plt.show()
    
    return(actual_control_normed)

def simd_sanger_heatmap_byexp(actual_controls, actual_results, initial_bits, expected_bits, exp_names):
    # Show results arranged together by experiment; similar to the NGS version in 20200507
    # Updated 8/19/20 to account for bad single direction reads
    # Updated 9/17/20 to add a column for controls too
    
    # Process actual data    
    actual_control_fwd, actual_control_rev = actual_controls
    actual_experiment_fwd, actual_experiment_rev = actual_results
    
    # Controls
    def normalize_bits(results):
        results_normed = []
        for sample in results:
            sample_bits = []
            for x,y in sample:
                if x == "?" or y == "?":
                    sample_bits.append((x, y))
                else:
                    sample_bits.append((x/(x+y), y/(x+y)))
            results_normed.append(sample_bits)
        return(results_normed)
                
    actual_control_fwd_normed = normalize_bits(actual_control_fwd)
    actual_control_rev_normed = normalize_bits(actual_control_rev)
    actual_experiment_fwd_normed = normalize_bits(actual_experiment_fwd)
    actual_experiment_rev_normed = normalize_bits(actual_experiment_rev)
    
    # Too messy to do by list comp
    def combined_normalized_bits(fwd_normed, rev_normed):
        actual_normed = []
        for sample_fwd, sample_rev in zip(fwd_normed, rev_normed):
            sample_normed_fwd = []
            sample_normed_rev = []

            for fwd, rev in zip(sample_fwd, sample_rev):
                # Save values based on which (Fwd or Rev) is viable
                # Odd rows will be fwd, even will be rev
                if "?" not in fwd and "?" not in rev:
                    sample_normed_fwd.append((fwd[1]+rev[1])/2)
                    sample_normed_rev.append((fwd[1]+rev[1])/2)
                elif "?" not in fwd and "?" in rev:
                    sample_normed_fwd.append(fwd[1])
                    sample_normed_rev.append(-1)
                elif "?" in fwd and "?" not in rev:
                    sample_normed_fwd.append(-1)
                    sample_normed_rev.append(rev[1])
                else: 
                    sample_normed_fwd.append(-1)
                    sample_normed_rev.append(-1)

            actual_normed.append(sample_normed_fwd)
            actual_normed.append(sample_normed_rev)
        return(actual_normed)
    
    actual_control_normed = combined_normalized_bits(actual_control_fwd_normed, actual_control_rev_normed)
    actual_experiment_normed = combined_normalized_bits(actual_experiment_fwd_normed, actual_experiment_rev_normed)
    
    # Organize comparison data
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [initial_values, actual_control_normed, expected_values, actual_experiment_normed]
    all_labels = ["Initial", "Control", "Expected +1", "Experiment"]
    
    # Plot
    n_samples = len(initial_bits)
    plt.close()
    fig, axes = plt.subplots(1, 4, figsize=(15, n_samples))
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_under(color="gray", alpha=1.0)
    for i in [0,2]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = axes[i]

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
    for i in [1,3]: # Dealing with the actual data
        data = all_data[i]
        label = all_labels[i]
        ax = axes[i]

        handle = ax.imshow(data, aspect=0.5, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False) # 4 = number of total bits
        ax.set_xticklabels([])
        
        ax.set_yticks([x*2-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples*2), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
    # Turn on label for just first plot
    axes[0].set_yticklabels([x + " " + y for x,y in zip(exp_names, initial_bits)], minor=True)
    
    # Label Fwd and Rev reads for actual data
    axes[-1].yaxis.tick_right()
    axes[-1].set_yticklabels([["F", "R"][i%2] for i in range(n_samples*2)], minor=True)
    axes[-1].yaxis.set_label_position("right")
    axes[-1].set_ylabel("Read direction", labelpad=20, rotation=270)
    
    # Make heatmap legend
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
    cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
    cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.show()
    
    return

def simd_sanger_heatmap_initial_actual(actual_results, expected_bits, \
                                       exp_names, legend=True, plot_adjust=0.76):
    # Simple Sanger version, just for presentation purposes
    # Only the expected and the experimental result are plotted
    
    # Process actual data    
    actual_experiment_fwd, actual_experiment_rev = actual_results
    
    # Controls
    def normalize_bits(results):
        results_normed = []
        for sample in results:
            sample_bits = []
            for x,y in sample:
                if x == "?" or y == "?":
                    sample_bits.append((x, y))
                else:
                    sample_bits.append((x/(x+y), y/(x+y)))
            results_normed.append(sample_bits)
        return(results_normed)
                
    actual_experiment_fwd_normed = normalize_bits(actual_experiment_fwd)
    actual_experiment_rev_normed = normalize_bits(actual_experiment_rev)
    
    # Too messy to do by list comp
    def combined_normalized_bits(fwd_normed, rev_normed):
        actual_normed = []
        for sample_fwd, sample_rev in zip(fwd_normed, rev_normed):
            sample_normed_fwd = []
            sample_normed_rev = []

            for fwd, rev in zip(sample_fwd, sample_rev):
                # Save values based on which (Fwd or Rev) is viable
                # Odd rows will be fwd, even will be rev
                if "?" not in fwd and "?" not in rev:
                    sample_normed_fwd.append((fwd[1]+rev[1])/2)
                    sample_normed_rev.append((fwd[1]+rev[1])/2)
                elif "?" not in fwd and "?" in rev:
                    sample_normed_fwd.append(fwd[1])
                    sample_normed_rev.append(-1)
                elif "?" in fwd and "?" not in rev:
                    sample_normed_fwd.append(-1)
                    sample_normed_rev.append(rev[1])
                else: 
                    sample_normed_fwd.append(-1)
                    sample_normed_rev.append(-1)

            actual_normed.append(sample_normed_fwd)
            actual_normed.append(sample_normed_rev)
        return(actual_normed)
    
    actual_experiment_normed = combined_normalized_bits(actual_experiment_fwd_normed, actual_experiment_rev_normed)
    
    # Organize comparison data
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [expected_values, actual_experiment_normed]
    all_labels = ["Expected", "Actual"]
    
    # Plot
    n_samples = len(expected_bits)
    plt.close()
    fig, axes = plt.subplots(1, 2, figsize=(9, n_samples))
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_under(color="gray", alpha=1.0)
    i = 0
    # We'll deal with the unknown bits later
    data = all_data[i]
    label = all_labels[i]
    ax = axes[i]

    handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
    ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
    ax.set_xticklabels([])

    ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
    ax.set_yticks(range(n_samples), minor=True)
    ax.set_yticklabels([])

    ax.grid(which="major", color="k", linestyle='-', linewidth=2)
    ax.set_title(label, y=1.02)
        
    i = 1 # Dealing with the actual data
    data = all_data[i]
    label = all_labels[i]
    ax = axes[i]

    handle = ax.imshow(data, aspect=0.5, cmap=current_cmap, vmin=0, vmax=1)
    ax.set_xticks([x + 0.5 for x in range(3)], minor=False) # 4 = number of total bits
    ax.set_xticklabels([])

    ax.set_yticks([x*2-0.5 for x in range(n_samples)], minor=False)
    ax.set_yticks(range(n_samples*2), minor=True)
    ax.set_yticklabels([])

    ax.grid(which="major", color="k", linestyle='-', linewidth=2)
    ax.set_title(label, y=1.02)

    # Turn on label for just first plot
    axes[0].set_yticklabels([x + " " + y for x,y in zip(exp_names, expected_bits)], minor=True)
    
    # Label Fwd and Rev reads for actual data
    axes[-1].yaxis.tick_right()
    axes[-1].set_yticklabels([["F", "R"][i%2] for i in range(n_samples*2)], minor=True)
    axes[-1].yaxis.set_label_position("right")
    axes[-1].set_ylabel("Read direction", labelpad=20, rotation=270)
    
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(right=plot_adjust)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
        cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.show()
    
    return

def cell_by_cell_barplot(listOfReadBits, title=None, legend=True):
    # Makes the Sanger-like plots for value by each cell
    # Data must be presorted, just a simple list of (readname, bit)
    colors = ["#f8971f", "#00a9b7", "gray"]
    possibleValues = ["0", "1", "?"]
    xPos=np.arange(4)
    width = 0.25
    
    # Count up bits for each cell and normalize
    nReads = len(listOfReadBits)
    relativeDigitCount = [{x: 0 for x in possibleValues} for _ in range(4)]
    for readname, bits in listOfReadBits:
        for digitindex, value in enumerate(bits):
            relativeDigitCount[digitindex][value] += 1
    relativeDigitCountNormalized = [{key: value/nReads for key, value in digit.items()} \
                                    for digit in relativeDigitCount]
    
    # Plot the bits
    plt.close()
    plt.figure(figsize=(8,5))
    for i, key in enumerate(possibleValues):
        yVals = [x[key] for x in relativeDigitCountNormalized]
        plt.bar(xPos+i*width, yVals, width, color=colors[i], label=key)
        for j, y in enumerate(yVals):
            plt.text(j+i*width, y+0.025, key, ha="center", color=colors[i])
    if legend:
        plt.legend(bbox_to_anchor=(1,0.5), loc=6, edgecolor="black", title="Value")
    plt.xticks(xPos+width, ["Cell {}".format(x) for x in range(1, 5)])
    plt.ylabel("Normalized read count")
    plt.xlabel("Digits")
    if title != None:
        plt.title(title, y=1.025)
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.show()
    
    return

def cell_by_cell_barplot_multiple(listOfListsOfReadBits, titles=None, \
                                  layout_shape=(2,2), suptitle=None, title_y=1.02):
    # Makes the Sanger-like plots for value by each cell
    # Data must be presorted, just a simple list of (readname, bit)
    # Puts all figures together as subplots into the specified shape
    # Goes down the column first before across
    
    # layout_shape elements must be >1 for indices to work
    if 1 in layout_shape or 0 in layout_shape:
        print("Can only take layout shapes of more than 1 row and column. Try again.")
        return
    
    colors = ["#f8971f", "#00a9b7", "gray"]
    possibleValues = ["0", "1", "?"]
    xPos=np.arange(4)
    width = 0.25
    nHeight, nWidth = layout_shape
    
    plt.close()
    fig, axes = plt.subplots(nHeight, nWidth, sharey=True, figsize=(8*nWidth, 6*nHeight))
    
    for k, (listOfReadBits,title) in enumerate(zip(listOfListsOfReadBits, titles)): 
        axX, axY = k%nHeight, math.floor(k/nHeight)
        
        # Count up bits for each cell and normalize
        nReads = len(listOfReadBits)
        relativeDigitCount = [{x: 0 for x in possibleValues} for _ in range(4)]
        for readname, bits in listOfReadBits:
            for digitindex, value in enumerate(bits):
                relativeDigitCount[digitindex][value] += 1
        relativeDigitCountNormalized = [{key: value/nReads for key, value in digit.items()} \
                                        for digit in relativeDigitCount]
            
        # Plot the bits
        for i, key in enumerate(possibleValues):
            yVals = [x[key] for x in relativeDigitCountNormalized]
            axes[axX][axY].bar(xPos+i*width, yVals, width, color=colors[i], label=key)
            for j, y in enumerate(yVals):
                axes[axX][axY].text(j+i*width, y+0.025, key, ha="center", color=colors[i])
        axes[axX][axY].set_xticks(xPos+width)
        axes[axX][axY].set_xticklabels(["Cell {}".format(x) for x in range(1, 5)])
        if axY == 0:
            axes[axX][axY].set_ylabel("Normalized read count")
        axes[axX][axY].set_xlabel("Digits")
        if titles != None:
            axes[axX][axY].set_title(title + "\nn = {}".format(nReads), y=1.025)
        axes[axX][axY].set_ylim(0, 1.1)
        axes[axX][axY].grid(True)
      
    customLegend = [Patch(facecolor=c, label=str(x)) \
                    for c,x in zip(colors, possibleValues)]

    axes[0][nHeight-1].legend(handles=customLegend, loc=1,\
              edgecolor="black", title="Value")

    if suptitle != None:
        plt.suptitle(suptitle, y=title_y, fontsize=24)
    fig.tight_layout()
    plt.show()
    
    return


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


def SIMD_barcode_distribution_inset(dataListofTuples, bitsofInterest, \
                                    readstoBarcodesandResults, title, \
                                    expectedN, nSIMDBarcodes=15, ncolor=0,\
                                   correspondingBitvals=None):
    # A plot of a subset of the data from one sample, to investigate 
    # reads with the wrong bits
    # Can be used like an "inset" plot, zooming into one subset
    # Not including bit values with ? (undetermined) digits in them
    # Includes reads with only barcode determined
    
    readcountsofInterestbySIMD = {x:0 for x in range(1, nSIMDBarcodes+1)}
    for readname, bitvalue in dataListofTuples:
        if bitvalue == bitsofInterest:
            info = readstoBarcodesandResults[readname]
            simdf, simdr = info["SIMD"]
            if simdf == simdr: 
                readcountsofInterestbySIMD[simdf] += 1
            elif simdf == 0 or simdr == 0:
                readcountsofInterestbySIMD[max(simdf, simdr)] += 1
            else: # Neither are determined; don't add
                continue
                
    data = [readcountsofInterestbySIMD[x] for x in range(1, nSIMDBarcodes+1)]
    normedData = [100*x/sum(data) for x in data]
    xPos = range(1, nSIMDBarcodes+1)
    width=0.75
    colors=["#f8971f", "#ffd600", "#a6cd57", "#579d42", "#00a9b7"]
    plt.close()
    plt.figure(figsize=(8,5))
    plt.bar(xPos, normedData, width, color=colors[ncolor])
    plt.title("{} ({} reads with bits {})".format(title, \
                sum(data), bitsofInterest),\
                 fontsize=16, y=1.05)
    plt.ylabel("% of reads in subset ".format(bitsofInterest))
    if correspondingBitvals == None:
        plt.xticks(xPos, xPos)
        plt.xlabel("SIMD Barcode")
    else: # Format to include 
        xtickswithbits = ["{}\n{}".format(x,y) for x,y in zip(xPos, correspondingBitvals)]
        plt.xticks(xPos, xtickswithbits, fontsize=10)
        plt.xlabel("SIMD barcode")
    plt.ylim(0,100)
    plt.gca().get_xticklabels()[xPos.index(expectedN)].set_color('red')
    plt.gca().get_xticklabels()[xPos.index(expectedN)].set_weight("bold")
    plt.grid(True)
    plt.show()
    return

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

def simd_ngs_heatmap_rule110(barcode_to_determined_bits, sample_index_to_barcodes,\
                             titles, initial_bits, expected_bits, legend=True, plot_adjust=0.76,\
                             mainTitle=None):
    # Simplest version of the viridis heatmap, just takes expected and 
    # the experimental outcome
    # Note that these plots take only both-NGS barcode reads
    # For presentation purposes. Can toggle legend on/off
    # The # of ? bits are indicated by size of rectangle
    # New 7/21/21: Makes n plots: initial, exp1, exp2, exp3
    # And takes their specific titles too
    # sample_index_to_barcodes is a list of list of barcodes
    
    possible_values = ["0", "1", "?"]
    n_samples = len(expected_bits)
    n_experiments = len(sample_index_to_barcodes)
    
    # Set up variables for later plotting
    overall_actual_digit_count_normed = []
    overall_unknown_count_normed = []
    overall_reads_counted = []
    
    for sample_index_to_barcode in sample_index_to_barcodes:
        actual_digit_count_normed = []
        unknown_count_normed = []
        reads_counted = []

        for i, barcode_pair in enumerate(sample_index_to_barcode):
            fwd_id, rev_id = barcode_pair

            # Organize data into list of dicts
            n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
            reads_counted.append(n_reads)
            if n_reads == 0:
                barcode_digit_count_normed = [np.nan for _ in range(4)]
                barcode_unknown_count_normed = [0 for _ in range(4)]
            else:
                rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
                for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                    readname, bit = readname_bit[0], readname_bit[1]
                    for digit_index, value in enumerate(bit):
                        rel_digit_count[digit_index][value] += 1
                norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
                rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]

                barcode_digit_count_normed = []
                for v in rel_digit_count_normed:
                    try:
                        barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                    except ZeroDivisionError:
                        barcode_digit_count_normed.append(np.nan)
                barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]

            actual_digit_count_normed.append(barcode_digit_count_normed)
            unknown_count_normed.append(barcode_unknown_count_normed)
            
        overall_actual_digit_count_normed.append(actual_digit_count_normed)
        overall_unknown_count_normed.append(unknown_count_normed)
        overall_reads_counted.append(reads_counted)
        
    # Prepare initial and expected values
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [initial_values, expected_values] + overall_actual_digit_count_normed
    all_labels = ["Initial","Expected"] + titles
    n_data = len(all_labels)
    
    # Plot
    widthPerCtl = 0.95/(2+1.2*n_experiments)
    widthPerData = widthPerCtl*1.2
    plt.close()
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    fig = plt.figure(figsize=(4*n_data+1, n_samples))
    
    gs_ctrls = GridSpec(nrows=1, ncols=2)
    gs_ctrls.update(left=0.0, right=2*widthPerCtl, wspace=0)
    
    #Deal with initial and exp bits first
    for i in range(n_data)[:2]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_ctrls[i])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        if i == 0: ax.set_yticklabels(initial_bits, minor=True)
    
    # Real data
    gs_expts = GridSpec(nrows=1, ncols=n_experiments)
    gs_expts.update(left=2*widthPerCtl, right=2*widthPerCtl+n_experiments*widthPerData, wspace=0)
    
    for i in range(n_data)[2:]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_expts[i-2])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        # Add right hand labels to Actual as read numbers
        reads_counted = overall_reads_counted[i-2]
        ax.set_yticklabels(["n = {}".format(x) for x in reads_counted], minor=True)
        ax.yaxis.tick_right()

        # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction
        unknown_count_normed = overall_unknown_count_normed[i-2]
        for r, row in enumerate(unknown_count_normed):
            for c, value in enumerate(row):
                if value > 0:
                    half_edge = math.sqrt(value)/2
                    rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                    ax.add_patch(rect)
    
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(right=plot_adjust)
        cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
        cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    if mainTitle != None:
        plt.suptitle(mainTitle,y=1.02, fontsize=24)
        
    plt.show()
    
    return

def simd_ngs_heatmap_multi(barcode_to_determined_bits, sample_index_to_barcodes,\
                             titles, expected_bits, legend=True, plot_adjust=0.76,\
                             mainTitle=None):
    # Takes expected value, plus n plots
    
    possible_values = ["0", "1", "?"]
    n_samples = len(expected_bits)
    n_experiments = len(sample_index_to_barcodes)
    
    # Set up variables for later plotting
    overall_actual_digit_count_normed = []
    overall_unknown_count_normed = []
    overall_reads_counted = []
    
    for sample_index_to_barcode in sample_index_to_barcodes:
        actual_digit_count_normed = []
        unknown_count_normed = []
        reads_counted = []

        for i, barcode_pair in enumerate(sample_index_to_barcode):
            fwd_id, rev_id = barcode_pair

            # Organize data into list of dicts
            n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
            reads_counted.append(n_reads)
            if n_reads == 0:
                barcode_digit_count_normed = [np.nan for _ in range(4)]
                barcode_unknown_count_normed = [0 for _ in range(4)]
            else:
                rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
                for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                    readname, bit = readname_bit[0], readname_bit[1]
                    for digit_index, value in enumerate(bit):
                        rel_digit_count[digit_index][value] += 1
                norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
                rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]

                barcode_digit_count_normed = []
                for v in rel_digit_count_normed:
                    try:
                        barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                    except ZeroDivisionError:
                        barcode_digit_count_normed.append(np.nan)
                barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]

            actual_digit_count_normed.append(barcode_digit_count_normed)
            unknown_count_normed.append(barcode_unknown_count_normed)
            
        overall_actual_digit_count_normed.append(actual_digit_count_normed)
        overall_unknown_count_normed.append(unknown_count_normed)
        overall_reads_counted.append(reads_counted)
        
    # Prepare initial and expected values
    expected_values = [[int(x) for x in entry] for entry in expected_bits]
    all_data = [expected_values] + overall_actual_digit_count_normed
    all_labels = ["Expected"] + titles
    n_data = len(all_labels)
    
    # Plot
    widthPerCtl = 0.95/(1+1.4*n_experiments)
    widthPerData = widthPerCtl*1.4
    plt.close()
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    fig = plt.figure(figsize=(4*n_data+1, n_samples))
    
    gs_ctrls = GridSpec(nrows=1, ncols=1)
    gs_ctrls.update(left=0.0, right=widthPerCtl, wspace=0)
    
    #Deal with expected bits first
    for i in range(n_data)[:1]: # Keeps as list
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_ctrls[i])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        if i == 0: ax.set_yticklabels(expected_bits, minor=True)
    
    # Real data
    gs_expts = GridSpec(nrows=1, ncols=n_experiments)
    gs_expts.update(left=1*widthPerCtl, right=1*widthPerCtl+n_experiments*widthPerData, wspace=0)
    
    for i in range(n_data)[1:]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_expts[i-1])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        # Add right hand labels to Actual as read numbers
        reads_counted = overall_reads_counted[i-1]
        ax.set_yticklabels(["n={}".format(x) for x in reads_counted], minor=True, fontsize=12)
        ax.yaxis.tick_right()

        # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction
        unknown_count_normed = overall_unknown_count_normed[i-1]
        for r, row in enumerate(unknown_count_normed):
            for c, value in enumerate(row):
                if value > 0:
                    half_edge = math.sqrt(value)/2
                    rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                    ax.add_patch(rect)
    
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(right=plot_adjust)
        cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
        cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    if mainTitle != None:
        plt.suptitle(mainTitle,y=0.96, fontsize=24)
        
    plt.show()
    
    return

def simd_ngs_heatmap_vertical(barcode_to_determined_bits, sample_index_to_barcodes,\
                              initial_bits, row_labels, \
                              legend=True, plot_adjust=0.05):
    # Plots heatmaps vertically. Dimensions are encoded in
    # sample_index_to_barcodes - list of lists
    # Number of lists = ncols
    # Number of items in lists inside lists = nrows
    # initial_bits length = ncols
    
    possible_values = ["0", "1", "?"]
    nCols = len(sample_index_to_barcodes)
    nRows = len(sample_index_to_barcodes[0])
    
    # Set up variables for later plotting; same dimensions as input
    overall_actual_digit_count_normed = []
    overall_unknown_count_normed = []
    overall_reads_counted = []
    
    for barcode_this_column in sample_index_to_barcodes:
        actual_digit_count_normed = []
        unknown_count_normed = []
        reads_counted = []
        
        for i, barcode_pair in enumerate(barcode_this_column):
            fwd_id, rev_id = barcode_pair
    
            # Organize data into list of dicts
            n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
            reads_counted.append(n_reads)
            if n_reads == 0:
                barcode_digit_count_normed = [np.nan for _ in range(4)]
                barcode_unknown_count_normed = [0 for _ in range(4)]
            else:
                rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
                for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                    readname, bit = readname_bit[0], readname_bit[1]
                    for digit_index, value in enumerate(bit):
                        rel_digit_count[digit_index][value] += 1
                norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
                rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]

                barcode_digit_count_normed = []
                for v in rel_digit_count_normed:
                    try:
                        barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                    except ZeroDivisionError:
                        barcode_digit_count_normed.append(np.nan)
                barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]

            actual_digit_count_normed.append(barcode_digit_count_normed)
            unknown_count_normed.append(barcode_unknown_count_normed)
            
        overall_actual_digit_count_normed.append(actual_digit_count_normed)
        overall_unknown_count_normed.append(unknown_count_normed)
        overall_reads_counted.append(reads_counted)
        
    # Prepare initial and expected values
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    all_data = []
    for initial_value, normed_counts in \
            zip(initial_values, overall_actual_digit_count_normed):
        all_data.append([initial_value] + normed_counts)
    
    # Plot!
    plt.close()
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    fig = plt.figure(figsize=(4*nCols+1, nRows*2))
    colWidth = 0.95/(1.4*nCols)
    plotWidth = 0.95/(nCols)
    colGap = colWidth - plotWidth
    
    for i, colData in enumerate(all_data):
        gs_col = GridSpec(nrows=nRows+1, ncols=1) # initial vals
        gs_col.update(left=i*plotWidth, right=i*plotWidth+colWidth, \
                      wspace=0)
        
        for j, rowData in enumerate(colData):
            ax = fig.add_subplot(gs_col[j])
            handle = ax.imshow([rowData], cmap=current_cmap, vmin=0, vmax=1)
            
            ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
            ax.set_xticklabels([])

            ax.set_yticks([-0.5,0.5], minor=False) # always be 1 
            ax.set_yticks([0], minor=True)
            ax.set_yticklabels([])

            ax.grid(which="major", color="k", linestyle='-', linewidth=2)
            
            # Row label
            if i == 0: # only do for the leftmost row
                if j == 0: # Control
                    ax.set_yticklabels(["Initial\nValue"], minor=True)
                else:
                    ax.set_yticklabels([row_labels[j-1]], minor=True)
            
            # Title
            if j == 0: # only do for the top row
                ax.set_title(initial_bits[i], y=1.02)
            else: # Read count; best done with text than twinx
                reads_counted = overall_reads_counted[i][j-1]
                ax.text(3.5*1.04, 0, "n={}".format(reads_counted), \
                        ha="left" , va="center", fontsize=12)
                
                # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction
                unknown_count_normed = overall_unknown_count_normed[i][j-1]
                for c, value in enumerate(unknown_count_normed):
                    if value > 0:
                        half_edge = math.sqrt(value)/2
                        rect = matplotlib.patches.Rectangle(\
                                (c-half_edge,0-half_edge), 2*half_edge, 2*half_edge, color="white")
                        ax.add_patch(rect)
            
            # Arrow
            if j != nRows: # final row
                ax.set_xlabel("$\downarrow$", fontsize=20)
            
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(bottom=plot_adjust)
        cbar_ax = fig.add_axes([0, 0-0.025, 0.9, 0.025])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1],\
                                  orientation="horizontal")
        cbar_handle.ax.set_xticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(-0.02, 0, "Bit Value", ha="right", va="center")
    
    plt.show()
    return

def simd_ngs_heatmap_rule110_controls(barcode_to_determined_bits, sample_index_to_barcodes,\
                             titles, initial_bits, legend=True, plot_adjust=0.76,\
                             mainTitle=None):
    # Simplest version of the viridis heatmap, just takes expected and 
    # the experimental outcome
    # Note that these plots take only both-NGS barcode reads
    # For presentation purposes. Can toggle legend on/off
    # The # of ? bits are indicated by size of rectangle
    # New 7/21/21: Makes n plots: initial, exp1, exp2, exp3
    # And takes their specific titles too
    # sample_index_to_barcodes is a list of list of barcodes
    
    possible_values = ["0", "1", "?"]
    n_samples = len(initial_bits)
    n_experiments = len(sample_index_to_barcodes)
    
    # Set up variables for later plotting
    overall_actual_digit_count_normed = []
    overall_unknown_count_normed = []
    overall_reads_counted = []
    
    for sample_index_to_barcode in sample_index_to_barcodes:
        actual_digit_count_normed = []
        unknown_count_normed = []
        reads_counted = []

        for i, barcode_pair in enumerate(sample_index_to_barcode):
            fwd_id, rev_id = barcode_pair

            # Organize data into list of dicts
            n_reads = len(barcode_to_determined_bits[fwd_id][rev_id])
            reads_counted.append(n_reads)
            if n_reads == 0:
                barcode_digit_count_normed = [np.nan for _ in range(4)]
                barcode_unknown_count_normed = [0 for _ in range(4)]
            else:
                rel_digit_count = [{x: 0 for x in possible_values} for _ in range(4)]
                for readname_bit in barcode_to_determined_bits[fwd_id][rev_id]:
                    readname, bit = readname_bit[0], readname_bit[1]
                    for digit_index, value in enumerate(bit):
                        rel_digit_count[digit_index][value] += 1
                norm_factor = len(barcode_to_determined_bits[fwd_id][rev_id])
                rel_digit_count_normed = [{key: value/norm_factor for key, value in digit.items()} for digit in rel_digit_count]

                barcode_digit_count_normed = []
                for v in rel_digit_count_normed:
                    try:
                        barcode_digit_count_normed.append(v["1"]/(v["1"]+v["0"]))
                    except ZeroDivisionError:
                        barcode_digit_count_normed.append(np.nan)
                barcode_unknown_count_normed = [v["?"] for v in rel_digit_count_normed]

            actual_digit_count_normed.append(barcode_digit_count_normed)
            unknown_count_normed.append(barcode_unknown_count_normed)
            
        overall_actual_digit_count_normed.append(actual_digit_count_normed)
        overall_unknown_count_normed.append(unknown_count_normed)
        overall_reads_counted.append(reads_counted)
        
    # Prepare initial and expected values
    initial_values = [[int(x) for x in entry] for entry in initial_bits]
    all_data = [initial_values] + overall_actual_digit_count_normed
    all_labels = ["Initial"] + titles
    n_data = len(all_labels)
    
    # Plot
    widthPerCtl = 0.95/(1+1.2*n_experiments)
    widthPerData = widthPerCtl*1.2
    plt.close()
    current_cmap = matplotlib.cm.get_cmap("viridis")
    current_cmap.set_bad(color="gray")
    fig = plt.figure(figsize=(4*n_data+1, n_samples))
    
    gs_ctrls = GridSpec(nrows=1, ncols=1)
    gs_ctrls.update(left=0.0, right=1*widthPerCtl, wspace=0)
    
    #Deal with initial and exp bits first
    for i in range(n_data)[:1]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_ctrls[i])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        if i == 0: ax.set_yticklabels(initial_bits, minor=True)
    
    # Real data
    gs_expts = GridSpec(nrows=1, ncols=n_experiments)
    gs_expts.update(left=1*widthPerCtl, right=1*widthPerCtl+n_experiments*widthPerData, wspace=0)
    
    for i in range(n_data)[1:]:
        # We'll deal with the unknown bits later
        data = all_data[i]
        label = all_labels[i]
        ax = fig.add_subplot(gs_expts[i-1])

        handle = ax.imshow(data, cmap=current_cmap, vmin=0, vmax=1)
        ax.set_xticks([x + 0.5 for x in range(3)], minor=False)
        ax.set_xticklabels([])
        
        ax.set_yticks([x-0.5 for x in range(n_samples)], minor=False)
        ax.set_yticks(range(n_samples), minor=True)
        ax.set_yticklabels([])
        
        ax.grid(which="major", color="k", linestyle='-', linewidth=2)
        ax.set_title(label, y=1.02)
        
        # Add right hand labels to Actual as read numbers
        reads_counted = overall_reads_counted[i-1]
        ax.set_yticklabels(["n = {}".format(x) for x in reads_counted], minor=True)
        ax.yaxis.tick_right()

        # Draw contrasting rectangle in bits w/ unknown reads, area proportional to unknown fraction
        unknown_count_normed = overall_unknown_count_normed[i-1]
        for r, row in enumerate(unknown_count_normed):
            for c, value in enumerate(row):
                if value > 0:
                    half_edge = math.sqrt(value)/2
                    rect = matplotlib.patches.Rectangle((c-half_edge,r-half_edge), 2*half_edge, 2*half_edge, color="white")
                    ax.add_patch(rect)
    
    # Make heatmap legend
    if legend:
        fig.subplots_adjust(right=plot_adjust)
        cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
        cbar_handle = fig.colorbar(handle, cax=cbar_ax, ticks=[0, 0.5, 1])
        cbar_handle.ax.set_yticklabels([0, "50/50\nmix", "1"])
        cbar_handle.ax.text(0, 1.05, "Bit Value", ha="center", va="center")
    
    if mainTitle != None:
        plt.suptitle(mainTitle,y=1.02, fontsize=24)
        
    plt.show()
    
    return

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

def expected_error_by_value_old(barcode_to_determined_bits, data_indices, expected_bits, title,\
                            toss_undetermined=True, n_samples=8, n_bits=4):
    # Plots distribution of bit values (correlated 4-bit values) by error, reports expected error
    
    plt.close()
    plt.figure(figsize=(n_samples*2, 5))
    colors = ["#f8971f", "#ffd600", "#a6cd57", "#579d42", "#00a9b7"]
    labels = ["{} errors".format(k) for k in np.arange(n_bits+1)]
    label_x_pos = np.arange((n_bits+2)/2,n_samples*(n_bits+2), (n_bits+2))-1
    
    all_expected_error = []
    
    for i in range(n_samples):
        correct_value = expected_bits[i]
        fwd_idx, rev_idx = data_indices[i]
        bits_list = [x[1] for x in barcode_to_determined_bits[fwd_idx][rev_idx]]
        sorted_bits = {x:[] for x in range(n_bits+1)}
        for bit in bits_list:
            if toss_undetermined and "?" in bit:
                continue
            else:
                error = hamming_dist(bit, correct_value)
                sorted_bits[error] += bit
    
        data = [len(sorted_bits[x]) for x in range(n_bits+1)]
        total_n = sum(data)
        if total_n == 0: 
            normed_data = [0 for _ in range(n_bits+1)] 
            ### skip this one if no data; not a good long term solution though!
        else:
            normed_data = [x/total_n for x in data]
        
        # Plot data
        width = 0.75
        for j in range(n_bits+1):
            x = i*(n_bits+2)+j
            plt.bar(x, normed_data[j], width, color=colors[j%5], label=labels[j])
            
        expected_error = np.sum([n*x for n,x in enumerate(normed_data)])
        plt.text(i*(n_bits+2), max(normed_data)*1.05, round(expected_error,2), \
                 color='black', ha="center")
        all_expected_error.append(expected_error)
    
    legend_handles = []
    for j in range(n_bits+1):
        legend_handles.append(Patch(color=colors[j])) 
    plt.legend(legend_handles, labels, fontsize=14, loc=1, ncol=n_bits+1)
    
    
    plt.ylim(0,1.2)
    plt.title("{} Bit sequence error (Hamming distance) by expected value".format(title), y=1.05)
    plt.ylabel("Proportion of reads (Normalized)")
    plt.xlabel("Expected bit value")
    plt.xticks(label_x_pos, expected_bits, ha="center")
    plt.grid(True)
    plt.show()
        
    return(all_expected_error)

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

def matrix_heatmap(plottingMatrix, referenceMatrix, title=None, title_adjust=1.1, \
                   select_rows=None, select_cols=None, outlineCells=True, \
                   stablepoints=None, cycles=None, threshold=None):
    # Given a (probability) matrix, makes a heatmap out of it
    # referenceMatrix is like a mask on top of the plottingMatrix
    # to specify which cells are labeled with the value
    # outlineCells option draws a box around the "correct" values for emphasis
    # stablepoints and cycles are bit values for special coloring to indicate
    # which values are stable and which are part of a cycle
    # cycles should be list of lists (each list is a cycle)
    InitVals = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',\
                 '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    
    if select_rows != None:
        xInitVals = [InitVals[x] for x in select_rows]
    else:
        xInitVals = InitVals
        
    if select_cols != None:
        yInitVals = [InitVals[x] for x in select_cols]
    else:
        yInitVals = InitVals
        
    dimensions = np.shape(plottingMatrix)
    plt.close()
    fig = plt.figure(figsize=(dimensions[1]*10/16,dimensions[0]*10/16))
    ax = fig.add_subplot(111)
    
    ### For a more fair visual comparison
    vmaxValue = 1
    if True:
        posMax = []
        for vec in plottingMatrix:
            posMax.append(max(vec))
        maxValue = max(posMax)
        if maxValue <= 1: vmaxValue = 1
        else: vmaxValue = maxValue
    
    handle = ax.imshow(plottingMatrix, cmap="magma", vmin=0, vmax=vmaxValue) 
    ax.set_xticks(np.arange(dimensions[1]))
    ax.set_xticklabels(yInitVals, rotation=90,fontsize=12)
    ax.set_yticks(np.arange(dimensions[0]))
    ax.set_yticklabels(xInitVals,fontsize=12)
    ax.tick_params(axis='x', which='major', \
                    labelbottom = False, bottom=False, top = False, labeltop=True)
    ax.tick_params(axis="y", which="major", length=0)
    if stablepoints != None:
        for point in stablepoints:
            ax.get_xticklabels()[xInitVals.index(point)].set_color("#f8971f")
            ax.get_xticklabels()[xInitVals.index(point)].set_weight("bold")
            ax.get_yticklabels()[yInitVals.index(point)].set_color("#f8971f")
            ax.get_yticklabels()[yInitVals.index(point)].set_weight("bold")
    if cycles != None:
        for cycle in cycles:
            for point in cycle:
                ax.get_xticklabels()[xInitVals.index(point)].set_color("#00a9b7")
                ax.get_xticklabels()[xInitVals.index(point)].set_weight("bold")
                ax.get_yticklabels()[yInitVals.index(point)].set_color("#00a9b7")
                ax.get_yticklabels()[yInitVals.index(point)].set_weight("bold")
            
    # Add text
    for i in range(dimensions[1]):
        for j in range(dimensions[0]):
            value = plottingMatrix[j][i]
            if referenceMatrix[j][i] != 1 and threshold == None:
                continue
            elif referenceMatrix[j][i] == 1 or value >= threshold:
                if value > 0.7: textcolor = "k"
                else: textcolor = "w"
                text = ax.text(i, j, "{0:.2f}".format(value),
                           ha="center", va="center", color=textcolor, fontsize=11)
                if outlineCells and referenceMatrix[j][i] == 1: # expected answer
                    ax.add_patch(Rectangle((i+0.475, j-0.475), -0.95, 0.95, fill=False, edgecolor='white', lw=3.))
                    ax.add_patch(Rectangle((i+0.43, j-0.41), -0.85, 0.85, fill=False, edgecolor='black', lw=1.5))
    if cycles !=None:
        labelspace = 30
    else: labelspace = 5
    ax.set_ylabel("Input",labelpad=labelspace)
    ax.set_xlabel("Output",labelpad=labelspace)
    ax.xaxis.set_label_position("top")
    if title != None:
        plt.title(title, y=title_adjust, va="bottom", ha="center")
        
    cax = fig.add_axes([0.945, 0.15, 0.025, 0.7]) # Colorbar
    cbar = fig.colorbar(handle, cax=cax)
    cbar.ax.text(-0.25, 0.5, "Fraction of reads", ha="right", va="center", rotation=90)
    plt.show()
    
    return