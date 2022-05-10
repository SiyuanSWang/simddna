from simddna.NGSanalysis import *

if __name__ == "__main__":

	##### ===== Parameters specific to run ===== #####

	# Algorithm "rule110" or "counting"?
	Algorithm = "counting" 
	
	# List all expected registers in the data
	RegisterNames = ["Reg8"]

	# Names of paired sequencing files (extension must be .fastq)
	Filenames = ["SIMDCountR8_S2_L001_R1_001.fastq",\
				 "SIMDCountR8_S2_L001_R2_001.fastq"]

	# NGS barcodes used in this run, 1-indexed
	BarcodesN = [9,10]
	
	# Does the read region cover the entire register?
	FullCoverage = False

	# Number of threads
	NThreads = 16

	##### ===== End of run parameters ===== #####

	# Filter viable reads 
	AllViableReads = filter_viable_reads(Filenames, Algorithm, RegisterNames, NThreads=NThreads)

	# Iterate through each register
	for register in RegisterNames:
		print("Processing Register {}...".format(register))

		# Identify NGS/Sample barcodes
		NGSIdentifiedReads = organize_by_ngs_barcodes(BarcodesN, AllViableReads[register], \
			register, FullCoverage, NThreads=NThreads)

		# Identify SIMD/Initial value barcodes
		SIMDIdentifiedReads = organize_by_simd_barcodes(register, len(BarcodesN), AllViableReads[register],\
				NGSIdentifiedReads, FullCoverage, ChunkSize=10000, NThreads=NThreads)

		# Read out bit values
		BarcodeToDeterminedBits = collect_bit_values(Algorithm, register, AllViableReads[register], \
				SIMDIdentifiedReads, len(BarcodesN), FullCoverage, NThreads=NThreads)

	print("Finished processing {} for {} algorithm!".format(\
			" ".join(RegisterNames), Algorithm.capitalize()))

