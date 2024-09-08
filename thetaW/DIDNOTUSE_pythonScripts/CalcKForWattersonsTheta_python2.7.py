import sys
import os
import gzip

if len (sys.argv) !=3:
	print("Usage: python CalcKForWattersonsTheta_python2.7.py <input gzipped vcf file> <output seg site, missing, mono counts file>")
	exit(1)

#Define Arguments
gvcfFile = sys.argv[1]
SegSiteoutFilename = sys.argv[2]

#Define Variables
endOfLine = '\n'
tab = '\t'
segValues = ['1/1','0/1','1/0', '1|1','0|1','1|0'] #GT of SNPS
MonoRef = ['0/0', '0|0'] #GT of SNPS
Missing = ['./.', '.|.'] #GT of SNPS

#Open files
datafile = gzip.open(gvcfFile, 'rt') #read as text
outFile = open(SegSiteoutFilename, 'w')
outFile.write('#CHROM\tPOS\tINDV\n')

for line in datafile:
	if line.startswith("##"): #skip lines that start with # character
		pass
	elif line.startswith("#CHROM"):
		header = line.strip().split('\t')
	else:
		dataCol = line.strip().split('\t')
		if dataCol[3] !='N': #skip uncallable stuff

		#Parse the Format Column for Genotype and grab position column 	
			FCOLS=dataCol[8].split(':')
			FCOLSlen=len(FCOLS)	

		#Find the GT info and count up the number of fixed REF or ALT alleles
			for i in range(9, len(dataCol)): #Genotype Cols start at 9
				if len(dataCol[i].split(':')) == FCOLSlen:	
					if dataCol[i].split(':')[FCOLS.index("GT")] in MonoRef:	
						pass
					elif dataCol[i].split(':')[FCOLS.index("GT")] in Missing:	
						pass
					elif dataCol[i].split(':')[FCOLS.index("GT")] in segValues:	
						outFile.write(str(dataCol[0]) + tab + str(dataCol[1]) + tab + str(header[i]) + endOfLine) #only write out segregating sites
					else:
						print(str(dataCol[0]) + tab + str(dataCol[1]) + tab + str(header[i]) + tab + dataCol[i].split(':')[FCOLS.index("GT")])
				else:
					pass
				
		else:
			pass

outFile.close()
