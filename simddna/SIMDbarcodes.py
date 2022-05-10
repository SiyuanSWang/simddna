# These are the SIMD/random access displacement barcodes
# used for all registers
# Barcodes are the same across the two algorithms
# but vary between M13 register sub-addresses and the IDT register
# Each initial value has its own SIMD barcode
# These variables get loaded as a module in NGSanalysis.py

InitVals = ["0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",\
			"1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"]

SIMDBarcodes = {}

# All registers used in NGS experiments
AllRegisters = ["IDT", "Reg1", "Reg2", "Reg3", "Reg7","Reg8","Reg9"] 

# IDT Barcodes; 30 nts long
SIMDBarcodes["IDT"] = {}
SIMDBarcodes["IDT"]["0000"] = "TGAGGTTTTGTTTTGTGAGTGGAATGGATG"
SIMDBarcodes["IDT"]["0001"] = "TGTTGTATGTTAGTGTGAATTGATATGTTG"
SIMDBarcodes["IDT"]["0010"] = "TGGTGTTGATTGATGTGGTTATGGATATTG"
SIMDBarcodes["IDT"]["0011"] = "TGGGTTTGAATAATGTGGTAATGAGTTGTG"

SIMDBarcodes["IDT"]["0100"] = "TGGATTTATGAAGTGTGGAATTATAGTGTG"
SIMDBarcodes["IDT"]["0101"] = "TGAAATGTAGGAATGTGTTTAGTATGGTTG"
SIMDBarcodes["IDT"]["0110"] = "TGTAAAGGATGAGTGTGGGATGAAGATTTG"
SIMDBarcodes["IDT"]["0111"] = "TGGTTTTGTTATGTGTGTATGTTTGAGTTG"

SIMDBarcodes["IDT"]["1000"] = "TGTGGATTTGTATTGTGGGTGTGTATAGTG"
SIMDBarcodes["IDT"]["1001"] = "TGATAGTGAAAAGTGTGTAATGGTTTGTTG"
SIMDBarcodes["IDT"]["1010"] = "TGTGAGTTGGGTTTGTGTGGGTTTTGGGTG"
SIMDBarcodes["IDT"]["1011"] = "TGAGATTTTAGGGTGTGTGGTTGTTATGTG"

SIMDBarcodes["IDT"]["1100"] = "TGGGAGTTAAGGATGTGGTGGTTGGTAATG"
SIMDBarcodes["IDT"]["1101"] = "TGGATGTAGTTTGTGTGTAGTTGTGAAGTG"
SIMDBarcodes["IDT"]["1110"] = "TGATATTAGGAGTTGTGAGTAGGTATGATG"
SIMDBarcodes["IDT"]["1111"] = "TGAAGTTTGGTGGTGTGATTGAGTGAATTG"

# All M13 Registers; 15 nts long
for reg in AllRegisters[1:]:
	SIMDBarcodes[reg] = {} 
	SIMDBarcodes[reg]["0000"] = "CAAAACAAAACCTCA"
	SIMDBarcodes[reg]["0001"] = "CAACATATCAATTCA"
	SIMDBarcodes[reg]["0010"] = "CATTATTCAAACCCA"
	SIMDBarcodes[reg]["0011"] = "CACTTCATAAATCCA"

	SIMDBarcodes[reg]["0100"] = "CAACTCCTAATATCA"
	SIMDBarcodes[reg]["0101"] = "CAACCATACTAAACA"
	SIMDBarcodes[reg]["0110"] = "CATTCCTACATTTCA"
	SIMDBarcodes[reg]["0111"] = "CAAATCTTCATCCCA"

	SIMDBarcodes[reg]["1000"] = "CAAACACTCTATTCA"
	SIMDBarcodes[reg]["1001"] = "CAACTCAAACATACA"
	SIMDBarcodes[reg]["1010"] = "CATACCCTTTTCTCA"
	SIMDBarcodes[reg]["1011"] = "CATCATACCTACTCA"

	SIMDBarcodes[reg]["1100"] = "CAAAACTCTCTCTCA"
	SIMDBarcodes[reg]["1101"] = "CAAACCCAACTCACA"
	SIMDBarcodes[reg]["1110"] = "CATTCTCCCACCTCA"
	SIMDBarcodes[reg]["1111"] = "CATATCTAATCTCCA"

