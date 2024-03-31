import sys
import pysam
import os
import gzip
import argparse

def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="This script computes sliding window heterozygosity.")

    parser.add_argument(
        "--vcf", required=True,
        help="REQUIRED. input vcf file, this must have passed through filtering and QC so *NO* bad sites remain in the file")
    parser.add_argument(
        "--window_size", required=True,
        help="REQUIRED. window size ")
    parser.add_argument("--step_size", required=True,
                       help="REQUIRED. Name of output file.")
    parser.add_argument("--chromNum", required=True,
                       help="REQUIRED. chromosome of query VCF")
    parser.add_argument("--chromLengths", required=True,
                       help="REQUIRED. A two-column tab-delimited list of chromosomes (col 1) and their lengths (col 2) in the ref. genome")
    parser.add_argument("--outpath", required=True,
                       help="REQUIRED. output file path")
    args = parser.parse_args()
    return args

args = parse_args()

# open input file (gzipped VCF file), make sure the VCF file is indexed (if not, create index)

filename = args.vcf

VCF = gzip.open(filename, 'rt')  # Use 'rt' for text mode in Python 3

if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# set remaining variables
window_size = int(args.window_size)
step_size = int(args.step_size)
QueryRegion = args.chromNum
outfilepath = args.outpath
het = ['0/1', '1/0', '0|1', '1|0']

# grab cromosome sizes
genomeSize = open(args.chromLengths, 'r')
chrom_size = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in genomeSize}
genomeSize.close()


# get list of samples
samples = []
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]:
            samples.append(i)
        break

# get first and last positions in chromosome
for line in VCF:
    if line[0] != '#':
        start_pos = int(line.strip().split()[1])
        end_pos = int(chrom_size[QueryRegion])
        break

# create output file
output = open(outfilepath + os.path.basename(filename) + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('QueryRegion\twindow_start\tsites_total\tsites_passing\tcalls_%s\tmissing_%s\thets_%s\thomRef_%s\thomAlt_%s\n' % ('\tcalls_'.join(samples), '\tmissing_'.join(samples), '\thets_'.join(samples), '\thomRef_'.join(samples), '\thomAlt_'.join(samples)) )

# Fetch a region, ignore sites that fail filters, tally total calls and heterozygotes
def snp_cal(QueryRegion, window_start, window_end):

    with gzip.open(filename, 'rt') as VCF:  # Re-open VCF within the function
        rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (QueryRegion, window_start, window_end), parser=pysam.asTuple()))

        sites_total, sites_passing = 0, 0
        calls = [0] * len(samples)
        hets = [0] * len(samples)
        missing = [0] * len(samples)
        homRef = [0] * len(samples)
        homAlt = [0] * len(samples)

        for line in rows:
            sites_total += 1
            sites_passing += 1
            for i in range(0, len(samples)):
                GT = line[i + 9]
                if GT[:1] != '.':
                    calls[i] += 1
                if GT[:3] in het:
                    hets[i] += 1

        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (QueryRegion,window_start,sites_total,sites_passing,'\t'.join(map(str,calls)),'\t'.join(map(str,missing)),'\t'.join(map(str,hets)),'\t'.join(map(str,homRef)),'\t'.join(map(str,homAlt))) )


# initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# calculate stats for window, update window start and end positions, repeat to end of chromosome
while window_end <= end_pos:

    if window_end < end_pos:
        snp_cal(QueryRegion,window_start,window_end)

        window_start = window_start + step_size
        window_end = window_start + window_size - 1

    else:
        snp_cal(QueryRegion,window_start,window_end)
        break

else:
    window_end = end_pos
    snp_cal(QueryRegion,window_start,window_end)


# close files and exit
VCF.close()
output.close()

exit()
