import sys
import os


if len (sys.argv) !=3:
	print "Usage: python create_SFS.py <input vcf file> <output seg site, missing, mono counts file>"
	exit(1)

vcfFilename = sys.argv[1]
SegSiteoutFilename = sys.argv[2]

#Counts
#COUNT_INVCF_LINES=0
#COUNT_PRINT=0

#Lists and Dictionaries
segValues = ['1/1','0/1','1/0'] #check to see that REF and ALT positions are valid alleles

if(os.path.isfile(vcfFilename)):

	with open(vcfFilename, 'r') as dataFile:
		with open(SegSiteoutFilename, 'w') as outFile:
			outFile.write('Chrom\tPos\tSegSite\n')

			for line in dataFile:
				if line.startswith("#"): #skip lines that begin with a hash
					pass
				else:
					SegSiteCount=0.0
					line = line.strip().split('\t')	

					#Get Derived Allele count per site
					for dataCol in line[9:]:
						if (dataCol in segValues):
							SegSiteCount=1
						elif dataCol == '0/0':
							pass
						elif dataCol == './.':
							pass
						else:
							print '['+dataCol+']', line[1]
							sys.exit('BAD GENOTYPE')
					#print line[1], SegSiteCount 
					outFile.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(SegSiteCount) + '\n')	

