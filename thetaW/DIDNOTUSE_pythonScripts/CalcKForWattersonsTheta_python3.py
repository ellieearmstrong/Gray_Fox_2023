import sys
import gzip

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 CalcKForWattersonsTheta_python3.py <input gzipped vcf file> <output seg site, missing, mono counts file>")
        sys.exit(1)

    # Define Arguments
    gvcfFile = sys.argv[1]
    SegSiteoutFilename = sys.argv[2]

    # Define Variables
    endOfLine = '\n'
    tab = '\t'
    segValues = ['1/1', '0/1', '1/0', '1|1', '0|1', '1|0']  # GT of SNPs
    MonoRef = ['0/0', '0|0']  # GT of SNPs
    Missing = ['./.', '.|.']  # GT of SNPs

    # Open files
    with gzip.open(gvcfFile, 'rt') as datafile, open(SegSiteoutFilename, 'w') as outFile:
        outFile.write('#CHROM\tPOS\tINDV\n')

        for line in datafile:
            if line.startswith("##"):  # skip lines that start with # character
                continue
            elif line.startswith("#CHROM"):
                header = line.strip().split('\t')
            else:
                dataCol = line.strip().split('\t')
                if dataCol[3] != 'N':  # skip uncallable stuff

                    # Parse the Format Column for Genotype and grab position column
                    FCOLS = dataCol[8].split(':')
                    FCOLSlen = len(FCOLS)

                    # Find the GT info and count up the number of fixed REF or ALT alleles
                    for i in range(9, len(dataCol)):  # Genotype Cols start at 9
                        if len(dataCol[i].split(':')) == FCOLSlen:
                            genotype = dataCol[i].split(':')[FCOLS.index("GT")]
                            if genotype in MonoRef:
                                continue
                            elif genotype in Missing:
                                continue
                            elif genotype in segValues:
                                outFile.write(f"{dataCol[0]}{tab}{dataCol[1]}{tab}{header[i]}{endOfLine}")  # only write out segregating sites
                            else:
                                print(f"{dataCol[0]}{tab}{dataCol[1]}{tab}{header[i]}{tab}{genotype}")
                        else:
                            continue

if __name__ == "__main__":
    main()
