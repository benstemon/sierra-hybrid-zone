import argparse
import os
import gzip


def skipHeader(line):
    """ For a line in VCF file, if line is not a header line, extract columns into list """
    if isinstance(line, bytes):
        line = line.decode('utf-8')  # Decode bytes into string using utf-8 encoding
    cols = line.rstrip('\n').split('\t')  # Remove newline character and split by tab
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

def main(args):
    info_file = open(args.info_file, 'r')
    vcf = gzip.open(args.vcf, 'rb')
    summaryfile = open(args.summaryfile, 'w')
    outvcf = gzip.open(args.outvcf, 'wt')
    headerfile = open(args.headerfile, 'r')
    summaryfile.write('chr\tbp\tafd\n')

    for hline in headerfile.readlines():
        outvcf.write(hline)

    allele_freq_diff = args.allele_freq_diff
    p1_mincount = args.p1count
    p2_mincount = args.p2count

    p1_indices = []
    p2_indices = []
    for line in info_file:
        cols = line.replace('\n','').split('\t')
        if cols[0] != 'sample':
            if cols[1] == '0':
                p2_indices.append(int(cols[2]))
            elif cols[1] == '1':
                p1_indices.append(int(cols[2]))

    for line in vcf:
        cols = skipHeader(line)
        if cols != 'header':
            LG = str(cols[0])
            bp = int(cols[1])
            p1_count = 0
            p2_count = 0
            p1_geno_counts = [0, 0, 0]
            p2_geno_counts = [0, 0, 0]
            category = ''
            for e in range(len(p1_indices)):
                geno = cols[9 + p1_indices[e]].split(':')[0]
                if geno != './.':
                    p1_count += 1
                    if geno == '0/0' or geno == '0|0':
                        p1_geno_counts[0] += 1
                    elif geno == '0/1' or geno == '0|1' or geno == '1|0':
                        p1_geno_counts[1] += 1
                    elif geno == '1/1' or geno == '1|1':
                        p1_geno_counts[2] += 1

            for l in range(len(p2_indices)):
                geno = cols[9 + p2_indices[l]].split(':')[0]
                if geno != './.':
                    p2_count += 1
                    if geno == '0/0' or geno == '0|0':
                        p2_geno_counts[0] += 1
                    elif geno == '0/1' or geno == '0|1' or geno == '1|0':
                        p2_geno_counts[1] += 1
                    elif geno == '1/1' or geno == '1|1':
                        p2_geno_counts[2] += 1

            if p1_count >= p1_mincount and p2_count >= p2_mincount:
                p1_ref_freq = (float(p1_geno_counts[0]) + 0.5*float(p1_geno_counts[1]))/float(p1_count)
                p2_ref_freq = (float(p2_geno_counts[0]) + 0.5*float(p2_geno_counts[1]))/float(p2_count)
                if abs(p1_ref_freq - p2_ref_freq) >= allele_freq_diff:
                    afd = p1_ref_freq - p2_ref_freq
                    summaryfile.write(str(LG) + '\t' + str(bp) + '\t' + str(afd) + '\n')
                    outvcf.write(line.decode('utf-8'))                
    info_file.close()
    summaryfile.close()
    outvcf.close()
    headerfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files.")
    parser.add_argument("-f", "--info_file", help="Path to the info file")
    parser.add_argument("-i", "--vcf", help="Path to the input VCF file")
    parser.add_argument("-s", "--summaryfile", help="Path to the output summary file")
    parser.add_argument("-o", "--outvcf", help="Path to the output VCF file")
    parser.add_argument("-m", "--headerfile", help="Path to the header file")
    parser.add_argument("-a", "--allele_freq_diff", type=float, default=1.0, help="Allele frequency difference threshold. Default is fixed differences ")
    parser.add_argument("--p1count", type=float, default=10, help="Minimum number of P1 samples without missing data needed to consider a site")
    parser.add_argument("--p2count", type=float, default=10, help="Minimum number of P2 samples without missing data needed to consider a site")
    args = parser.parse_args()
    main(args)

# info-file = tab-delimited file with three columns: sample name, syndrome, and index in vcf. Counting is base-0. "Syndrome" denotes p1 and p2, where p1 = "1" and p2 = "0".
