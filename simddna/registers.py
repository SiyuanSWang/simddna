# All constant regions and cells of registers

from simdlib.SIMDlib import dna_complement

# Lengths of registers (distance between barcodes)
NGSBuffer = 8 # Buffer around region to look for SIMD barcode
NGSBarcodeLength = 11
NGSBarcodeLocations = {}
DistNGSBarcodes = {"IDT":230, "Reg1":222, "Reg2":209, "Reg3":312, "Reg7":312, "Reg8":315, "Reg9":303}
for register in DistNGSBarcodes.keys():
	NGSBarcodeLocations[register] = (DistNGSBarcodes[register] + NGSBarcodeLength-NGSBuffer, \
									DistNGSBarcodes[register] + 2*NGSBarcodeLength+NGSBuffer)

SIMDBuffer = 5
SIMDBarcodeLengthIDT = 30
SIMDBarcodeLengthM13 = 15
SIMDBarcodeLocations = {}
SIMDFwdReadPositions = {"IDT":196, "Reg1":203, "Reg2":190, "Reg3":203, "Reg7":293, "Reg8":296, "Reg9":284}
for register in SIMDFwdReadPositions.keys():
	if register[:3] == "IDT": BarcodeLength = SIMDBarcodeLengthIDT
	else: BarcodeLength = SIMDBarcodeLengthM13
	SIMDBarcodeLocations[register] = {"F":(SIMDFwdReadPositions[register]-SIMDBuffer, \
											SIMDFwdReadPositions[register]+BarcodeLength+SIMDBuffer),\
									"R":(26-SIMDBuffer, 26+BarcodeLength+SIMDBuffer)} 
									# Always in the same position in the reverse read

# Constant regions are split up differently based on which algorithm is used
AllConstantRegionsFwd = {"counting":{}, "rule110":{}}

# Binary Counting Algorithm
AllConstantRegionsFwd["counting"]["IDT"] = ["TCGTCGGCAGCGTCGATTAAGGATTTGAAGTTTTG",\
											"GAAATGGAATGATGATTGTGAGATTGGTGTTTGG",\
											"TATGTAGTGGAGTAGAGAGAGTTAGGTTG",\
											"AGGAGAAAGGTTTGTTGATAGATGTTTGAGTAAA",\
											"GTGTGTGTAGGAGATGTAAAGGTATGTGTGAGTTGGAAGTGAGAGTTGT",
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["counting"]["Reg3"] = ["TCGTCGGCAGCGTCCAAACTACAACGCCTGTAGC",\
											"TTCCACAGACAGCCCTCATAGTTAGCGTAACGA",\
											"CTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAA",\
											"GAATTTTCTGTATGGGATTTTGCTAAACAACTTTC",\
											"ACAGTTTCAGCGGAGTGAGAATAGAAAGGAACAACTAAAGGAATTGCGAATA",
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["counting"]["Reg7"] = ["TCGTCGGCAGCGTCTTTAGGCAGAGGCATTTTCGAG", \
											"CAGTAATAAGAGAATATAAAGTACCGACAAAAGGTAAAGTAATTCTGTCCAGACG",\
											"CGACAATAAACAACATGTTCAGCTAATGCAGAACGCGCCTGTTTATCAACAAT",\
											"GATAAGTCCTGAACAAGAAAAATAATATCCCATCCTAATTTACGAGCATGTAGAAACC",\
											"ATCAATAATCGGCTGTCTTTCCTTATCATTCCAAGAACGGGTATTAAACCAAGTACCGCACTCATCGAGAACAAGC",\
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["counting"]["Reg8"] = ["TCGTCGGCAGCGTCTTGTTTGGATTATACTTCTGAA",\
											"AATGGAAGGGTTAGAACCTACCATATCAAAATTATTTGCACGTAAAACAGAAATAAAGAAA",\
											"TGCGTAGATTTTCAGGTTTAACGTCAGATGAATATACAGTAACAGTACCTTTTACA",\
											"CGGGAGAAACAATAACGGATTCGCCTGATTGCTTTGAATACCAAGTTACAAAA",\
											"CGCGCAGAGGCGAATTATTCATTTCAATTACCTGAGCAAAAGAAGATGATGAAACAAACATCAAGAAAACAAAAT",
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["counting"]["Reg9"] = ["TCGTCGGCAGCGTCACAATATTTTTGAATGGCTA",\
											"TAGTCTTTAATGCGCGAACTGATAGCCCTAAAACATCGCCATTAAAAATACCGAACGAACC",\
											"CCAGCAGAAGATAAAACAGAGGTGAGGCGGTCAGTATTAACACCGCCTGC",\
											"ACAGTGCCACGCTGAGAGCCAGCAGCAAATGAAAAATCTAAAGCATCACC",\
											"TGCTGAACCTCAAATATCAAACCCTCAATCAATATCTGGTCAGTTGGCAAATCAACAGTTGAAAGGAATTGAGG",\
											"CCGAGCCCACGAGAC"]											

# Rule 110 Algorithm
AllConstantRegionsFwd["rule110"]["IDT"] = ["GATTAAGGATTTGAAGTTTTGTGAAATGGAATGATGA",\
										"TGTGAGATTGGTGTTTGGTTATGTAGTGGAG",\
										"AGAGAGAGTTAGGTTGAAGGAGAAAGGTTTGTTGA",\
										"AGATGTTTGAGTAAATGTGTGTGTAGGAGATG", \
										"AAAGGTATGTGTGAGTTGGAAGTGAGAGTTGT",\
										"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["rule110"]["Reg1"] = ["AGAACCGGATATTCATTACCCAAATCAACGTAACAAAG",\
											"TGCTCATTCAGTGAATAAGGCTTGCCCTGACG",\
											"GAAACACCAGAACGAGTAGTAAATTGGGCTTGAG",\
											"TGGTTTAATTTCAACTTTAATCATTGTGAATTACCTT", \
											"TGCGATTTTAAGAACTGGCTCATTATACCAGTC",\
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["rule110"]["Reg2"] = ["GCTTGATACCGATAGTTGCGCCGACAATGACAACAACC", \
											"TCGCCCACGCATAACCGATATATTCGGTCG", \
											"TGAGGCTTGCAGGGAGTTAAAGGCCGCT", \
											"TTGCGGGATCGTCACCCTCAGCAGCGAAAGAC", \
											"GCATCGGAACGAGGGTAGCAACGGCTACAGAGG", \
											"CCGAGCCCACGAGAC"]

AllConstantRegionsFwd["rule110"]["Reg3"] = ["CAAACTACAACGCCTGTAGCATTCCACAGACAGCCCTCA", \
											"AGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTT", \
											"CAGACGTTAGTAAATGAATTTTCTGTATGGGA", \
											"TTTGCTAAACAACTTTCAACAGTTTCAGCGGAG", \
											"GAGAATAGAAAGGAACAACTAAAGGAATTGCGAATA",\
											"CCGAGCCCACGAGAC"]

# Cell sequences are also split by the algorithm
FwdCells = {"counting":{},"rule110":{}}
FwdBits = {"counting":{},"rule110":{}}

# Binary Counting Algorithm
FwdCells["counting"]["IDT"] = ["AGGAGAAAGGTTTGTTGATAGATGTTTGAGTAAATGTGTGTGTAGGAGATGTAAAGGTATGTG",\
							"TATGTAGTGGAGTAGAGAGAGTTAGGTTGAAGGAGAAAGGTTTGTTGATAGATGTTTGAGTAAA",\
							"GAAATGGAATGATGATTGTGAGATTGGTGTTTGGTTATGTAGTGGAGTAGAGAGAGTTAGGTTG",\
							"GATTAAGGATTTGAAGTTTTGTGAAATGGAATGATGATTGTGAGATTG"]
FwdBits["counting"]["IDT"] = ["GTAAA"+"T"+"GTGTG",\
							"GGTTG"+"A"+"AGGAG",\
							"TTTGG"+"T"+"TATGT",\
							"TTTTG"+"T"+"GAAAT"]

FwdCells["counting"]["Reg3"] = ["TTCTGTATGGGATTTTGCTAAACAACTTTCAACAGTTTCAGCGGAGTGAGAATAGAAAGGA",\
								"AGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGGGATTTTGCTAAACAA",\
								"CACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAG",\
								"CAAACTACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGA"]
FwdBits["counting"]["Reg3"] = ["CTTTC"+"A"+"ACAGT",\
							"GTAAA"+"T"+"GAATT",\
							"AACGA"+"T"+"CTAAA",\
							"GTAGC"+"A"+"TTCCA"]

FwdCells["counting"]["Reg7"] = ["TAGAAACCAATCAATAATCGGCTGTCTTTCCTTATCATTCCAAGAACGGGTATTAAACCA",\
								"TATCAACAATAGATAAGTCCTGAACAAGAAAAATAATATCCCATCCTAATTTACGAGCATG",\
								"GTAATTCTGTCCAGACGACGACAATAAACAACATGTTCAGCTAATGCAGAACGCGCCTGTT",\
								"TTTAGGCAGAGGCATTTTCGAGCCAGTAATAAGAGAATATAAAGTACCGACAAAAGGTAAA"]
FwdBits["counting"]["Reg7"] = ["AAACC"+"A"+"ATCAA",\
							"ACAAT"+"A"+"GATAA",\
							"AGACG"+"A"+"CGACA",\
							"TCGAG"+"C"+"CAGTA"]

FwdCells["counting"]["Reg8"] = ["CAAGTTACAAAATCGCGCAGAGGCGAATTATTCATTTCAATTACCTGAGCAAAAGAAGATGAT",\
								"AGTAACAGTACCTTTTACATCGGGAGAAACAATAACGGATTCGCCTGATTGCTTTGAATAC",\
								"CACGTAAAACAGAAATAAAGAAATTGCGTAGATTTTCAGGTTTAACGTCAGATGAATATAC",\
								"TTGTTTGGATTATACTTCTGAATAATGGAAGGGTTAGAACCTACCATATCAAAATTATTTG"]
FwdBits["counting"]["Reg8"] = ["CAAAA"+"T"+"CGCGC",\
							"TTACA"+"T"+"CGGGA",\
							"AGAAA"+"T"+"TGCGT",\
							"CTGAA"+"T"+"AATGG"]

FwdCells["counting"]["Reg9"] = ["AGCATCACCTTGCTGAACCTCAAATATCAAACCCTCAATCAATATCTGGTCAGTTGGCA",\
								"TATTAACACCGCCTGCAACAGTGCCACGCTGAGAGCCAGCAGCAAATGAAAAATCTAA",\
								"GCCATTAAAAATACCGAACGAACCACCAGCAGAAGATAAAACAGAGGTGAGGCGGTCAG",\
								"ACAATATTTTTGAATGGCTATTAGTCTTTAATGCGCGAACTGATAGCCCTAAAACATC"]
FwdBits["counting"]["Reg9"] = ["TCACC"+"T"+"TGCTG",\
							"CCTGC"+"A"+"ACAGT",\
							"GAACC"+"A"+"CCAGC",\
							"GGCTA"+"T"+"TAGTC"]

# Rule 110 Algorithm
FwdCells["rule110"]["IDT"] = ["AGATGTTTGAGTAAATGTGTGTGTAGGAGATGTAAAGGTATGTGTGAGTTGGAAGTGAGAGTTGT",\
							"AGAGAGAGTTAGGTTGAAGGAGAAAGGTTTGTTGATAGATGTTTGAGTAAATGTGTGTGTAGGAGATG",\
							"TGTGAGATTGGTGTTTGGTTATGTAGTGGAGTAGAGAGAGTTAGGTTGAAGGAGAAAGGTTTGTTGA",\
							"GATTAAGGATTTGAAGTTTTGTGAAATGGAATGATGATTGTGAGATTGGTGTTTGGTTATGTAGTGGAG"]
FwdBits["rule110"]["IDT"] = ["AGATG"+"T"+"AAAGG",\
							"GTTGA"+"T"+"AGATG",\
							"TGGAG"+"T"+"AGAGA",\
							"GATGA"+"T"+"TGTGA"]

FwdCells["rule110"]["Reg1"] = ["TGGTTTAATTTCAACTTTAATCATTGTGAATTACCTTATGCGATTTTAAGAACTGGCTCATTATACCAGTC",\
							"GAAACACCAGAACGAGTAGTAAATTGGGCTTGAGATGGTTTAATTTCAACTTTAATCATTGTGAATTACCTT",\
							"TGCTCATTCAGTGAATAAGGCTTGCCCTGACGAGAAACACCAGAACGAGTAGTAAATTGGGCTTGAG",\
							"AGAACCGGATATTCATTACCCAAATCAACGTAACAAAGCTGCTCATTCAGTGAATAAGGCTTGCCCTGACG"]
FwdBits["rule110"]["Reg1"] = ["ACCTT"+"A"+"TGCGA",\
							"TTGAG"+"A"+"TGGTT",\
							"TGACG"+"A"+"GAAAC",\
							"CAAAG"+"C"+"TGCTC"]

FwdCells["rule110"]["Reg2"] = ["TTGCGGGATCGTCACCCTCAGCAGCGAAAGACAGCATCGGAACGAGGGTAGCAACGGCTACAGAGG",\
							"TGAGGCTTGCAGGGAGTTAAAGGCCGCTTTTGCGGGATCGTCACCCTCAGCAGCGAAAGAC",\
							"TCGCCCACGCATAACCGATATATTCGGTCGCTGAGGCTTGCAGGGAGTTAAAGGCCGCT",\
							"GCTTGATACCGATAGTTGCGCCGACAATGACAACAACCATCGCCCACGCATAACCGATATATTCGGTCG"]
FwdBits["rule110"]["Reg2"] = ["AAGAC"+"A"+"GCATC",\
							"CCGCT"+"T"+"TTGCG",\
							"GGTCG"+"C"+"TGAGG",\
							"CAACC"+"A"+"TCGCC"]

FwdCells["rule110"]["Reg3"] = ["TTTGCTAAACAACTTTCAACAGTTTCAGCGGAGTGAGAATAGAAAGGAACAACTAAAGGAATTGCGAATA",\
							"CAGACGTTAGTAAATGAATTTTCTGTATGGGATTTTGCTAAACAACTTTCAACAGTTTCAGCGGAG",\
							"AGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGGGA",\
							"CAAACTACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTT"]
FwdBits["rule110"]["Reg3"] = ["CGGAG"+"T"+"GAGAA",\
							"TGGGA"+"T"+"TTTGC",\
							"TCTTT"+"C"+"CAGAC",\
							"CCTCA"+"T"+"AGTTA"]

# Prepare reverse cells and bits
RevCells = {"counting":{},"rule110":{}}
RevBits = {"counting":{},"rule110":{}}
for algorithm in FwdCells.keys():
	for register in FwdCells[algorithm].keys():
		RevCells[algorithm][register] = [dna_complement(x) for x in FwdCells[algorithm][register]]
		RevBits[algorithm][register] = [dna_complement(x) for x in FwdBits[algorithm][register]]
