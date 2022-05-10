import gzip
import shutil
import pickle
import math
from itertools import product
from multiprocessing import Pool
from datetime import datetime
import os
from Bio import SeqIO
from functools import partial
from simddna.SIMDlib import *
from simddna.SIMDbarcodes import * # SIMD Barcode sequences
from simddna.registers import * # Register-affiliated locations and sequences
from simddna.NGSbarcodes import * # NGS barcode sequences (11-2)

def make_constant_region_mutations(RegisterNames, AllConstantRegionsFwd, \
		NToleratedMM=1):
	AllConstantRegionsFwdMuts = {x: [] for x in RegisterNames}
	for register in RegisterNames:
		for domain in AllConstantRegionsFwd[register]:
			subset = set([domain])
			subset = subset.union(set(list(all_n_subs(domain,NToleratedMM))))
			subset = subset.union(set(list(all_n_subs(domain,NToleratedMM))))
			subset = subset.union(set(list(all_n_dels(domain,NToleratedMM))))
			subset = subset.union(set(all_n_ins(domain,NToleratedMM)))
			AllConstantRegionsFwdMuts[register].append(list(subset))
			print("Added all allowed mutants of forward domain" + \
					" {} for {}.".format(domain, register))
		print("Finished forward domains for {}\n".format(register))

	AllConstantRegionsRevMuts = {x: [] for x in RegisterNames}
	for register in RegisterNames:
		for domain in AllConstantRegionsFwd[register]:
			domainrev = dna_complement(domain)
			subset = set([domainrev])
			subset = subset.union(set(list(all_n_subs(domainrev,NToleratedMM))))
			subset = subset.union(set(list(all_n_subs(domainrev,NToleratedMM))))
			subset = subset.union(set(list(all_n_dels(domainrev,NToleratedMM))))
			subset = subset.union(set(all_n_ins(domainrev,NToleratedMM)))
			AllConstantRegionsRevMuts[register].append(list(subset))
			print("Added all allowed mutants of reverse domain" + \
					" {} for {}.".format(domain, register))
		print("Finished reverse domains for {}\n".format(register))
	return(AllConstantRegionsFwdMuts, AllConstantRegionsRevMuts)

def prepare_directories(RegisterNames, Filenames):
	for register in RegisterNames:
			if not os.path.exists(register):
				os.mkdir(register)

	# Uncompress .gz files
	for filename in Filenames:
		if not os.path.exists(filename):
			print("File {} doesn't exist. Extracting from gzip...".format(\
					filename))
			try: 
				with gzip.open(filename+".gz", 'rb') as f_in:
				    with open(filename, 'wb') as f_out:
				        shutil.copyfileobj(f_in, f_out)
			except FileNotFoundError:
				print("Missing the gzip file {} in this directory.".format(\
						filename+".gz"))
				os._exit(0)
	return()

def viable_read_check(rawReads, RegisterNames, AllConstantRegionsFwdMuts, \
		AllConstantRegionsRevMuts, NConsecutiveAcceptable=3):
	readID, PEID, seq, qscore = rawReads
	for register in RegisterNames:
		if PEID[0] == "1": # forward read
			constantRegionsMuts = AllConstantRegionsFwdMuts[register]
		else: # reverse read
			constantRegionsMuts = AllConstantRegionsRevMuts[register]

		passed=[False for _ in range(len(AllConstantRegionsFwdMuts[register]))] 
		newlyIdentifiedRegister = False

		for i, domainSet in enumerate(constantRegionsMuts):
			# Search through all the constant domains
			for domain in domainSet: # all the variants
				if domain in seq: # Found it
					passed[i] = True
					break # next domain

		if True in passed:
			newlyIdentifiedRegister = register

		# check that we have the correct # of consecutive "True"
		accepted = False
		count = 0
		for value in passed:
			if value == True:
				count += 1
			else: # chain broke; reset
				count = 0
			if count == NConsecutiveAcceptable: 
				accepted = True
				break

		if accepted: # break out of bigger loop
			break 

	if accepted:
		return(newlyIdentifiedRegister, readID, PEID, seq, qscore)
	else: # Maybe register is identified, but otherwise not a qualified read
		return((None, readID, PEID, seq, qscore))

def filter_viable_reads(Filenames, Algorithm, RegisterNames, NThreads=16, \
		saveChunkSize=1000000):
	regsWithFiles = []
	for register in RegisterNames:
		firstFileNamebyReg = "{}/viable-reads-0mil.pkl".format(register)
		regsWithFiles.append(os.path.exists(firstFileNamebyReg))

	if True in regsWithFiles: # Stop and report which are done
		print("At least 1 viable reads file found for:")
		for register, found in zip(RegisterNames, regsWithFiles):
			if found: print(register)
			# Ends process

		# Load viable reads
		nMil = 0 
		AllViableReads = {r:{} for r in RegisterNames}
		for register in RegisterNames:
			while True: 
				try:
					with open("{}/viable-reads-{}mil.pkl".format(\
							register, nMil), "rb") as f:
						tempAllViableReads = pickle.load(f)
						print("Loaded viable reads " + \
								"{} to {} mil ".format(nMil, nMil+1) + \
								"for register {}.".format(register))
					for read, readDict in tempAllViableReads.items():
						AllViableReads[register][read] = readDict
					print("Added reads from " + \
							"{} to {} mil to".format(nMil, nMil+1) + \
							"AllViableReads for register {}.".format(register))
					nMil += 1
				except FileNotFoundError:
					break
			print("Added {} reads for register {}.".format(\
					len(AllViableReads[register]), register))

	else: # All false; proceed with viable reads extraction
		print("Didn't find viable reads file. Creating one...\n")

		AllConstantRegionsFwdMuts, AllConstantRegionsRevMuts = \
			make_constant_region_mutations(RegisterNames, \
					AllConstantRegionsFwd[Algorithm])

		prepare_directories(RegisterNames, Filenames)

		viable_read_check_params = partial(viable_read_check, \
				RegisterNames=RegisterNames, \
				AllConstantRegionsFwdMuts=AllConstantRegionsFwdMuts, \
				AllConstantRegionsRevMuts=AllConstantRegionsRevMuts)

		AllViableReads = {r:{} for r in RegisterNames} 
		for file in Filenames: # Fastq, not gzip
			newReads = 0
			print("\nReading file {}...".format(file))

			rawReads = []
			for record in SeqIO.parse(file, "fastq"):
				readname = record.name # Note: no "@"
				sequence = str(record.seq)
				score = record.letter_annotations["phred_quality"]
				direction = record.description.split(" ")[1]

				rawReads.append((readname, direction, sequence, score))

			# Process in parallel
			with Pool(NThreads) as p:
				ViableResults = p.map(viable_read_check_params, rawReads)

			for result in ViableResults:
				if result[0] == None:
					continue # Don't include
				else:
					register, readID, PEID, seq, qscore = result
					try: 
						AllViableReads[register][readID][PEID] = (seq, qscore)
					except KeyError: # Is new read, or partner didn't make it
						AllViableReads[register][readID] = {}
						AllViableReads[register][readID][PEID] = (seq, qscore)
					newReads += 1

			print("Included {} reads out of {} from file {}.".format(\
					newReads, len(rawReads), file))

		# Second pass to catch all missed partners
		for file in Filenames:
			# Need to catch reads in second file too
			rescuedReads = 0
			print("Going back for a second pass to add " + \
					"paired reads from file {}".format(file))

			for record in SeqIO.parse(file, "fastq"):
				readname = record.name # Note: no "@"
				sequence = str(record.seq)
				score = record.letter_annotations["phred_quality"]
				direction = record.description.split(" ")[1]

				# Check if paired read is already in the dictionary
				for register in RegisterNames:
					if readname in AllViableReads[register]:
						if direction not in AllViableReads[register][readname]:
							AllViableReads[register][readname][direction] = \
									(sequence, score)
							rescuedReads += 1
							break
		
			print("Added {} reads from {} ".format(rescuedReads, file) + \
					"because their paired read passed.")

		for r in RegisterNames:
			print("Register {}: Total {} reads out of {} are viable!".format(\
					r, len(AllViableReads[r]), len(rawReads)))
		
		for register in RegisterNames:
			readnameList = list(AllViableReads[register].keys())
			for i in range(math.floor(\
					len(AllViableReads[register])/saveChunkSize)+1):
				subViableReads = {}
				for read in readnameList[saveChunkSize*i:saveChunkSize*(i+1)]:
					subViableReads[read] = AllViableReads[register][read]
				with open("{}/viable-reads-{}mil.pkl".format(register, i), \
						"wb") as f:
					pickle.dump(subViableReads, f)
				print("Saved viable reads for register {}".format(register) + \
						" from {} to {} million to file.".format(i,i+1))

	return(AllViableReads)	

def prepare_ngs_barcodes(BarcodesN, mutsAllowedNGS=2):

	# Organize barcode numbering relative to expected barcodes
	SelectedBarcodes = []
	for i in BarcodesN:
		SelectedBarcodes.append(NGSBarcodes[i-1])

	# Note: NGS barcodes will be using new indices!
	NGSOrderedFwdBarcodes = [] 
	for barcode in SelectedBarcodes:
		subset = set([barcode]) # Use set to avoid repeats
		subset = subset.union(set(list(all_n_subs(barcode,mutsAllowedNGS))))
		subset = subset.union(set(list(all_n_dels(barcode,mutsAllowedNGS))))
		subset = subset.union(set(list(all_n_ins(barcode,mutsAllowedNGS))))
		NGSOrderedFwdBarcodes.append(list(subset))

	NGSOrderedRevBarcodes = [] 
	for barcode in SelectedBarcodes:
		barcoderev = dna_complement(barcode)
		subset = set([barcoderev]) # Use set to avoid repeats
		subset = subset.union(set(list(all_n_subs(barcoderev,mutsAllowedNGS))))
		subset = subset.union(set(list(all_n_dels(barcoderev,mutsAllowedNGS))))
		subset = subset.union(set(list(all_n_ins(barcoderev,mutsAllowedNGS))))
		NGSOrderedRevBarcodes.append(list(subset))

	return(NGSOrderedFwdBarcodes, NGSOrderedRevBarcodes)

def determine_ngs_barcodes(readDicts, FullCoverage, DistNGSBarcodes, \
		NGSOrderedFwdBarcodes, NGSOrderedRevBarcodes):
	readname, pairedReadDict = readDicts
	if FullCoverage:
		barcodeCalls = NGS_barcode_identifier_full_coverage(pairedReadDict, \
			NGSOrderedFwdBarcodes, NGSOrderedRevBarcodes, \
			pred3BarcodeRange=DistNGSBarcodes)
	else:
		barcodeCalls = NGS_barcode_identifier_incomplete_coverage(\
			pairedReadDict, NGSOrderedFwdBarcodes, NGSOrderedRevBarcodes, \
			pred1stBarcodeRange=(0,13), barcodeLength=NGSBarcodeLength)
	fwdBarcode, revBarcode = barcodeCalls

	return(readname, fwdBarcode, revBarcode)

def organize_by_ngs_barcodes(BarcodesN, allViableReads, \
		RegisterName, FullCoverage, NThreads=16, NSplit=100):
	# allViableReads is a subset of viable reads specific to the register
	NGSIdentifiedReads = {x: {y: [] for y in range(len(BarcodesN) + 1)} \
								for x in range(len(BarcodesN) + 1)}

	try:
		for n,m in product(range(len(BarcodesN) + 1), repeat=2):
			with open("{}/ngs_barcode_calls_F{}_R{}.pkl".format(\
					RegisterName, n, m), "rb") as f:
				NGSIdentifiedReads[n][m] = pickle.load(f)
		print("Loaded barcode calls. Barcodes have already been identified.")

	except FileNotFoundError: 
		print("Could not find barcode calls. Decoding barcodes now...")

		NGSOrderedFwdBarcodes, NGSOrderedRevBarcodes = \
				prepare_ngs_barcodes(BarcodesN, mutsAllowedNGS=2)

		# Iterate for each register
		print("Starting NGS barcode identification...")

		readnamesList = list(allViableReads.keys())
		chunksize = math.floor(len(readnamesList)/NSplit) # Get the mod reads
		if chunksize == 0 or len(readnamesList) == 0:
			raise ValueError("No viable reads found. Check the ealier step?")

		determine_ngs_barcodes_params = partial(determine_ngs_barcodes, \
				FullCoverage=FullCoverage, \
				DistNGSBarcodes=NGSBarcodeLocations[RegisterName], \
				NGSOrderedFwdBarcodes=NGSOrderedFwdBarcodes, \
				NGSOrderedRevBarcodes=NGSOrderedRevBarcodes)

		for i in range(NSplit):
			chunk = readnamesList[i*chunksize:(i+1)*chunksize]
			if i > 0:
				print("Looked at {} out of {} reads so far".format(\
						i*chunksize, len(readnamesList)))

			readDicts = [(x,allViableReads[x]) for x in chunk]

			print("##### ===== Starting a new pool chunk! =====")
			with Pool(NThreads) as p: 
				barcodeResults = p.map(determine_ngs_barcodes_params, readDicts)

			for result in barcodeResults:
				# Don't care what the result is, save it just the same
				readname, fwdBarcode, revBarcode = result

				fwdBarcodeSave, revBarcodeSave = fwdBarcode, revBarcode
				if fwdBarcode == None:
					fwdBarcodeSave = 0
				if revBarcode == None:
					revBarcodeSave = 0

				NGSIdentifiedReads[fwdBarcodeSave][revBarcodeSave].append(\
						readname)

			# Print the last read, just to get a peep into decoding
			print("Identified NGS Barcodes " + \
					"Fwd {} and Rev {} for {} read {}.".format(\
						fwdBarcode, revBarcode, RegisterName, readname))
			print("Finished pool chunk at "+str(datetime.now())+"\n")

		# Save separate file for each barcode combo
		for n,m in product(range(len(BarcodesN) + 1), repeat=2):
			with open("{}/ngs_barcode_calls_F{}_R{}.pkl".format(\
					RegisterName, n, m), "wb") as f:
				pickle.dump(NGSIdentifiedReads[n][m], f)
				print("Saved barcode " + \
						"Fwd {} Rev {} for register {} to file.".format(\
							n,m,RegisterName))

		print("Wrote all barcode calls to file.")

	return(NGSIdentifiedReads)

def prepare_simd_barcodes(RegisterName, mutsAllowedSIMD=1):
	# List of lists of permitted mutated barcodes
	SIMDOrderedFwdBarcodes = []
	for bit in InitVals:
		barcode = SIMDBarcodes[RegisterName][bit]
		subset = set([barcode]) # to avoid repeats
		subset = subset.union(set(list(all_n_subs(barcode,mutsAllowedSIMD))))
		subset = subset.union(set(list(all_n_dels(barcode,mutsAllowedSIMD))))
		subset = subset.union(set(list(all_n_ins(barcode,mutsAllowedSIMD))))
		SIMDOrderedFwdBarcodes.append(list(subset))

	SIMDOrderedRevBarcodes = []
	for bit in InitVals:
		barcode = SIMDBarcodes[RegisterName][bit]
		barcoderev = dna_complement(barcode)
		subset = set([barcoderev]) # to avoid repeats
		subset = subset.union(set(list(all_n_subs(barcoderev,mutsAllowedSIMD))))
		subset = subset.union(set(list(all_n_dels(barcoderev,mutsAllowedSIMD))))
		subset = subset.union(set(list(all_n_ins(barcoderev,mutsAllowedSIMD))))
		SIMDOrderedRevBarcodes.append(list(subset))

	return(SIMDOrderedFwdBarcodes, SIMDOrderedRevBarcodes)

def determine_simd_barcodes(pairedReadDict, FullCoverage, \
		SIMDOrderedFwdBarcodes, SIMDOrderedRevBarcodes, RegisterName):

	if FullCoverage: # Read covers entire register, barcode is in both reads
		SIMDbarcode = SIMD_barcode_identifier_full_coverage(\
				pairedReadDict, SIMDOrderedFwdBarcodes, SIMDOrderedRevBarcodes,\
				predRangeinFwdRead=SIMDBarcodeLocations[RegisterName]["F"], \
				predRangeinRevRead=SIMDBarcodeLocations[RegisterName]["R"]) 
	
	else: # Barcode is only in the reverse read
		SIMDbarcode = SIMD_barcode_identifier_incomplete_coverage(\
			pairedReadDict, SIMDOrderedFwdBarcodes, SIMDOrderedRevBarcodes, \
			predRangeinRevRead=SIMDBarcodeLocations[RegisterName]["R"])
		# Just use the default barcode ranges (set for Reg8)
	
	return(SIMDbarcode) 

def organize_by_simd_barcodes(RegisterName, TotalBarcodesNGS, AllViableReads, \
		NGSIdentifiedReads, FullCoverage, ChunkSize=10000, NThreads=16, \
		TotalBarcodesSIMD=16):
	try: # Don't decode SIMD if it's already done
		with open("{}/simd_barcode_calls.pkl".format(RegisterName), "rb") as f:
			SIMDIdentifiedReads = pickle.load(f)
			print("SIMD barcodes have already been identified.")

	except FileNotFoundError:
		print("Couldn't find target barcode call file. Decoding SIMD barcodes.")	
		# Sort by identified SIMD; a is SIMD barcode by Fwd read, b is Rev
		SIMDIdentifiedReads = {x: {y: {} \
								for y in range(TotalBarcodesNGS + 1)} \
								for x in range(TotalBarcodesNGS + 1)}

		SIMDOrderedFwdBarcodes, SIMDOrderedRevBarcodes = \
				prepare_simd_barcodes(RegisterName)

		for n, m in product(range(TotalBarcodesNGS + 1), repeat=2):
			readsByNGSBarcode = NGSIdentifiedReads[n][m]
			printText = "Barcode {}-{}".format(n,m)

			if FullCoverage:
				readsWithSIMDBarcode = \
						{a: {b: [] for b in range(TotalBarcodesSIMD + 1)} \
							for a in range(TotalBarcodesSIMD + 1)} 
							# [(readname, (SIMD F read, SIMD R read))...]
			else:
				readsWithSIMDBarcode = \
						{a: [] for a in range(TotalBarcodesSIMD + 1)}
						# [(readname, (SIMD read))...]

			determine_simd_barcodes_params = partial(determine_simd_barcodes, \
					FullCoverage=FullCoverage,\
					SIMDOrderedFwdBarcodes=SIMDOrderedFwdBarcodes, \
					SIMDOrderedRevBarcodes=SIMDOrderedRevBarcodes, \
					RegisterName=RegisterName)

			NSplit = math.ceil(len(readsByNGSBarcode)/ChunkSize)
			for i in range(NSplit):
				readsByNGSBarcodeSplit = \
					readsByNGSBarcode[i * ChunkSize:(i+1) * ChunkSize]
				nReads = len(readsByNGSBarcodeSplit)
				readsByNGSBarcodeSplitDict = \
					[AllViableReads[read] for read in readsByNGSBarcodeSplit]

				print("##### ===== Starting a new pool chunk for {}, ".format(\
						printText) + "{} out of {} ===== ".format(i+1, NSplit))

				with Pool(NThreads) as p:
					SIMDbarcodeResults = p.map(\
							determine_simd_barcodes_params, \
							readsByNGSBarcodeSplitDict)

				for readname, simdN in zip(readsByNGSBarcodeSplit, \
						SIMDbarcodeResults):
					if FullCoverage:
						simdF, simdR = simdN
						readsWithSIMDBarcode[simdF][simdR].append(readname)
					else:
						readsWithSIMDBarcode[simdN].append(readname)

				print("Identified SIMD Barcode {} for read {}.".format(\
						simdN, readname))

				print("\tCompleted chunk {} out of {} at ".format(i+1,NSplit)+\
						str(datetime.now()) + "\n")

			SIMDIdentifiedReads[n][m] = readsWithSIMDBarcode

		# Don't need to save things by different NGS barcode - it's not too big
		with open("{}/simd_barcode_calls.pkl".format(RegisterName), "wb") as f:
			pickle.dump(SIMDIdentifiedReads, f)
			print("Wrote NGS and SIMD barcode calls to file.")

	return(SIMDIdentifiedReads)

def prepare_bit_sequences(Algorithm, RegisterName, upToNGaps=2):
	# Includes all versions of the bits sequence with up to N number of gaps
	FwdBitsGapped = []
	for fwdBit in FwdBits[Algorithm][RegisterName]:
		bitVars = []
		for n in range(upToNGaps+1):
			bitVars += n_noncontig_gaps(fwdBit, n)
		FwdBitsGapped.append(bitVars)

	RevBitsGapped = []
	for revBit in RevBits[Algorithm][RegisterName]:
		bitVars = []
		for n in range(upToNGaps+1):
			bitVars += n_noncontig_gaps(revBit, n)
		RevBitsGapped.append(bitVars)

	return(FwdBitsGapped, RevBitsGapped)

def determine_bit_values(seqDicts, FwdCells, RevCells, FwdBitsGapped, \
		RevBitsGapped, Algorithm):
	result = cell_bit_caller_phred(seqDicts, (FwdCells, RevCells), \
			(FwdBitsGapped, RevBitsGapped), algorithm=Algorithm)
	return("".join(result))

def collect_bit_values(Algorithm, RegisterName, AllViableReads, \
		SIMDIdentifiedReads, TotalBarcodesNGS, FullCoverage, \
		NThreads=16, TotalBarcodesSIMD=16):
	try: 
		BarcodeToDeterminedBits = {}
		for ngsF in range(TotalBarcodesNGS+1):
			with open("{}/interpreted-bits_ngsF-{}.pkl".format(\
					RegisterName,ngsF), "rb") as f:
				BarcodeToDeterminedBits[ngsF] = pickle.load(f)
				print("Loaded decoded bits for Register {} NGS Fwd {}".format(\
						RegisterName, ngsF))
		print("Bits have been identified for {}.".format(RegisterName))  

	except FileNotFoundError:
		print("Couldn't find bit-determined reads for " + \
				"{}.\nReading bits...".format(RegisterName))

		BarcodeToDeterminedBits = {ngsF:{ngsR:{simdF:{simdR:[] \
				for simdR in range(TotalBarcodesSIMD+1)} \
				for simdF in range(TotalBarcodesSIMD+1)} \
				for ngsR in range(TotalBarcodesNGS+1)} \
				for ngsF in range(TotalBarcodesNGS+1)}

		FwdBitsGapped, RevBitsGapped = \
				prepare_bit_sequences(Algorithm, RegisterName)

		determine_bit_values_params = partial(determine_bit_values, \
				FwdCells=FwdCells[Algorithm][RegisterName], \
				RevCells=RevCells[Algorithm][RegisterName], \
				FwdBitsGapped=FwdBitsGapped, RevBitsGapped=FwdBitsGapped, \
				Algorithm=Algorithm)

		# Parallelize bit decoding
		for ngsF in range(TotalBarcodesNGS+1):
			for ngsR in range(TotalBarcodesNGS+1):
				if FullCoverage: # two simd barcodes
					simdIterator = product(range(TotalBarcodesSIMD+1), repeat=2)
				else:
					simdIterator = range(TotalBarcodesSIMD+1)

				for simdN in simdIterator:
					if FullCoverage:
						simdF, simdR = simdN
						printText = "NGS {}-{}, SIMD {}-{}".format(\
								ngsF, ngsR, simdF, simdR)
						readnamesThisBarcode = \
								SIMDIdentifiedReads[ngsF][ngsR][simdF][simdR]
					else:
						printText = "NGS {}-{}, SIMD {}".format(\
								ngsF, ngsR, simdN)
						readnamesThisBarcode = \
								SIMDIdentifiedReads[ngsF][ngsR][simdN]

					if len(readnamesThisBarcode) == 0:
						continue

					nReads = len(readnamesThisBarcode)
					seqDicts = \
							[AllViableReads[readname] \
									for readname in readnamesThisBarcode]

					print("Decoding for {}, size {}, started at ".format(\
							printText, nReads) + str(datetime.now()))
					with Pool(NThreads) as p:
						resultBits = \
								p.map(determine_bit_values_params, seqDicts)

					packedReadnames = \
							[(x, y) for x,y in \
									zip(readnamesThisBarcode, resultBits)]
					
					peekRead, peekBits = packedReadnames[-1] # Sneak a peak
					print("Decoded read {} for {} with value {}".format(\
							peekRead, printText, peekBits))

					if FullCoverage:
						BarcodeToDeterminedBits[ngsF][ngsR][simdF][simdR] = \
								packedReadnames
					else:
						BarcodeToDeterminedBits[ngsF][ngsR][simdN] = \
								packedReadnames

			# Save it in parts by forward NGS barcode
			with open("{}/interpreted-bits_ngsF-{}.pkl".format(\
					RegisterName,ngsF),'wb') as f:
				pickle.dump(BarcodeToDeterminedBits[ngsF], f)
				print("Wrote decoded bits for " + \
						"NGS Fwd {} reads to file.".format(ngsF))

	return(BarcodeToDeterminedBits)

