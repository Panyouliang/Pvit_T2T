import sys
import argparse

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This is a script for extract top length isoform from annotation files
that including different isoforms infomation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2024-06-13, yyyy-mm-dd
======================================================================''')

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-paf1', metavar='.paf', type=str, required=True, help='Please input the PAF file')
parser.add_argument('-paf2', metavar='.paf', type=str, required=True, help='Please input the PAF file')
parser.add_argument('-genome', metavar='fasta', type=str, required=True, help='Please input the genome fasta file')
args = parser.parse_args()

#================================================================================


def readpaf(F,length):
    genome_range = {}
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split('\t')
            name = line[5]
            if int(line[10]) > length:
                if name not in genome_range:
                    genome_range[name] = [[int(line[7]),int(line[8])]]
                else:
                    genome_range[name].append([int(line[7]),int(line[8])])

    genome_range_merge = {}
    for k,v in genome_range.items():
        vsort = sorted(v, key=lambda x:int(x[1]))
        genome_range_merge[k] = merge_coordinates(vsort)

    return genome_range_merge

def merge_coordinates(coordinates):
    coordinates.sort(key=lambda x: x[0])
    merged = []
    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] <= current_end:
            current_end = max(current_end, coord[1])
        else:
            merged.append([current_start, current_end])
            current_start, current_end = coord
    merged.append([current_start, current_end])
    return merged


def stat_genome(F):
    statd = {}
    gseqd = {}
    chrs,length = '',0
    total_L = 0
    with open(F,'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if chrs != '':
                    statd[chrs] = length
                    gseqd[chrs] = ''.join(seq)
                chrs = line[1:].split()[0]
                length = 0
                seq = []
            else:
                length += len(line)
                total_L += len(line)
                seq.append(line)
        statd[chrs] = length
        gseqd[chrs] = ''.join(line)

    return statd, gseqd, total_L

def reverse_intervals(total_length, existing_intervals):
    existing_intervals.sort(key=lambda x: x[0])
    reversed_intervals = []
    current_start = 1
    for start, end in existing_intervals:
        if start > current_start:
            reversed_intervals.append([current_start, start - 1])
        current_start = end + 1

    if current_start <= total_length:
        reversed_intervals.append([current_start, total_length])
    return reversed_intervals


def target_range(target_loc, loc_list):
    target_start, target_end = target_loc
    target_length = target_end - target_start
    overlap_count = 0
    total_overlap_length = 0
    part_site = []

    for start,end in loc_list:
        if start <= target_start and end >= target_end:
            total_overlap_length = target_length
            part_site.append([target_start, target_end])
            break
        elif start < target_end and end > target_start:
            overlap_count += 1
            overlap_length = min(end, target_end) - max(start, target_start)
            total_overlap_length += overlap_length
            part_site.append([max(start, target_start), min(end, target_end)])

    overlap_ratio = total_overlap_length / target_length if target_length > 0 else 0

    return overlap_ratio, total_overlap_length, part_site




if __name__ == '__main__':

    Length = 0
    paf_range_minimap  = readpaf(args.paf1,Length) # minimap.paf
    paf_range_winnomap = readpaf(args.paf2,Length) # winnomap.paf
    len_stats, gseq, length = stat_genome(args.genome) # genome.fa

    miss_range_minimap = {}
    total_miss_len_mini = 0
    for chrs, loc_ in paf_range_minimap.items():
        reverse_loc_ = reverse_intervals(len_stats[chrs], loc_)
        miss_range_minimap[chrs] = reverse_loc_
        for site in reverse_loc_:
            total_miss_len_mini += site[1] - site[0]

    miss_range_winnomap = {}
    total_miss_len_wino = 0
    for chrs, loc_ in paf_range_winnomap.items():
        reverse_loc_ = reverse_intervals(len_stats[chrs], loc_)
        miss_range_winnomap[chrs] = reverse_loc_
        for site in reverse_loc_:
            total_miss_len_wino += site[1] - site[0]

    print(f"#Chrs\tMissName\tStart\tEnd")
    total_L = 0
    for chrs, loclist in miss_range_minimap.items():
        chrs_miss_L = 0
        num = 0
        for site in loclist:
            if chrs in miss_range_winnomap:
                ration, cut_len, part_site = target_range(site, miss_range_winnomap[chrs])
                total_L += cut_len
                chrs_miss_L += cut_len
                for loc in part_site:
                    if loc[1] - loc[0] > 80:
                        num += 1
                        print(f"{chrs}\t{chrs}-miss-{str(num)}\t{loc[0]}\t{loc[1]}")

