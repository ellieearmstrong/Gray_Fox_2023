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

# set variables
window_size = int(args.window_size)
step_size = int(args.step_size)
chromo = args.chromNum
outfilepath = args.outpath
het = ['0/1', '1/0', '0|1', '1|0']

#Chromosome sizes
#Chromosome size humans
#chromo_size={'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,'6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,'11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,'16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,'21':46709983,'22':50818468}

#Chromosome size dogs canfam3.1
#chromo_size={'chr1':122678785,'chr2':85426708,'chr3':91889043,'chr4':88276631,'chr5':88915250,'chr6':77573801,'chr7':80974532,'chr8':74330416,'chr9':61074082,'chr10':69331447,'chr11':74389097,'chr12':72498081,'chr13':63241923,'chr14':60966679,'chr15':64190966,'chr16':59632846,'chr17':64289059,'chr18':55844845,'chr19':53741614,'chr20':58134056,'chr21':50858623,'chr22':61439934,'chr23':52294480,'chr24':47698779,'chr25':51628933,'chr26':38964690,'chr27':45876710,'chr28':41182112,'chr29':41845238,'chr30':40214260,'chr31':39895921,'chr32':38810281,'chr33':31377067,'chr34':42124431,'chr35':26524999,'chr36':30810995,'chr37':30902991,'chr38':23914537,'chrX':123869142}

#Chromosome size dogs canfam4
#chromo_size={'chr1':123555469,'chr2':84977118,'chr3':92478159,'chr4':89533978,'chr5':89562246,'chr6':78111829,'chr7':81081196,'chr8':76405309,'chr9':61171409,'chr10':70642554,'chr11':74803998,'chr12':72969919,'chr13':64297965,'chr14':61111200,'chr15':64675283,'chr16':60359999,'chr17':65086665,'chr18':56472573,'chr19':55514801,'chr20':58626990,'chr21':51741255,'chr22':61573379,'chr23':53134597,'chr24':48565227,'chr25':51730345,'chr26':39256414,'chr27':46661788,'chr28':41732930,'chr29':42516734,'chr30':40643382,'chr31':39900754,'chr32':40224581,'chr33':32138516,'chr34':42397673,'chr35':28050405,'chr36':31223115,'chr37':30785715,'chr38':24802598,'chrX':124987930}

#Chromosome size gray foxes
chromo_size={'chr1':117572690,'chr2':102603151,'chr3':98526141,'chr4':96410777,'chr5':91265922,'chr6':90720172,'chr7':82759371,'chr8':79079208,'chr9':74450539,'chr10':73623094,'chr11':73040223,'chr12':72908785,'chr13':71581588,'chr14':66704989,'chr15':64818041,'chr16':63979375,'chr17':59774042,'chr18':64801261,'chr19':58983324,'chr20':63093800,'chr21':56962569,'chr22':54392966,'chr23':53374214,'chr24':52741366,'chr25':69993150,'chr26':49334718,'chr27':47458807,'chr28':62992706,'chr29':79160281,'chr30':41889812,'chr31':40660448,'chr32':92326701}

#Chromosome sizes arctic foxes
#chromo_size={'chr1':185856410,'chr2':179390407,'chr3':165895547,'chr4':152257076,'chr5':133749769,'chr6':134365958,'chr7':131533942,'chr8':134799312,'chr9':79894757,'chr10':116666219,'chr11':110415469,'chr12':85157011,'chr13':61872336,'chr14':36868461,'chr15':49723682,'chr16':61660327,'chr17':43367292,'chr18':47775060,'chr19':53523694,'chr20':38527645,'chr21':45731683,'chr22':31001485,'chr23':63747709,'chr24':58816956,'chrX':123034249}

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
        end_pos = chromo_size[chromo]
        break

# create output file
output = open(outfilepath + os.path.basename(filename) + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('chromo\twindow_start\tsites_total\tsites_passing\tcalls_%s\tmissing_%s\thets_%s\thomRef_%s\thomAlt_%s\n' % ('\tcalls_'.join(samples), '\tmissing_'.join(samples), '\thets_'.join(samples), '\thomRef_'.join(samples), '\thomAlt_'.join(samples)) )

# Fetch a region, ignore sites that fail filters, tally total calls and heterozygotes
def snp_cal(chromo, window_start, window_end):

    with gzip.open(filename, 'rt') as VCF:  # Re-open VCF within the function
        rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chromo, window_start, window_end), parser=pysam.asTuple()))

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

        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chromo,window_start,sites_total,sites_passing,'\t'.join(map(str,calls)),'\t'.join(map(str,missing)),'\t'.join(map(str,hets)),'\t'.join(map(str,homRef)),'\t'.join(map(str,homAlt))) )


# initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# calculate stats for window, update window start and end positions, repeat to end of chromosome
while window_end <= end_pos:

    if window_end < end_pos:
        snp_cal(chromo,window_start,window_end)

        window_start = window_start + step_size
        window_end = window_start + window_size - 1

    else:
        snp_cal(chromo,window_start,window_end)
        break

else:
    window_end = end_pos
    snp_cal(chromo,window_start,window_end)


# close files and exit
VCF.close()
output.close()

exit()
