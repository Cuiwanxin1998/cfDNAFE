import pysam
import pybedtools
from collections import defaultdict
import os
import itertools
import numpy as np
import pandas as pd
import math
from skmisc.loess import loess
import gzip
from bx.intervals.intersection import Intersecter, Interval
from functools import partial
import math
from collections import Iterable


class commonError(Exception):
    def __init__(self, message):
        self.message = message


def maxCore(nCore=None):
    if nCore > 16:
        return 16
        print("The thread number is forced to 16!")
    else:
        return nCore


def rmEndString(x, y):
    for item in y:
        if x.endswith(item):
            x = x.replace(item, "")

    return x


def isSoftClipped(cigar):
    """
    cigar information:
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    """
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False


def GCcontent(seq):
    nA = seq.count("a") + seq.count("A")
    nT = seq.count("t") + seq.count("T")
    nG = seq.count("g") + seq.count("G")
    nC = seq.count("c") + seq.count("C")
    percent_GC = (nG + nC) / (nA + nT + nG + nC)
    return percent_GC


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    reference:
        https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        # filter reads
        if read.is_unmapped or read.is_qcfail or read.is_duplicate:
            continue
        if not read.is_paired:
            continue
        if not read.is_proper_pair:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mate_is_unmapped:
            continue
        if read.rnext != read.tid:
            continue
        if read.template_length == 0:
            continue
        if isSoftClipped(read.cigar):
            continue

        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def reverse_seq(seq):
    r_seq = ''
    for i in seq:
        if i == 'A':
            r_seq = r_seq + 'T'
        elif i == 'T':
            r_seq = r_seq + 'A'
        elif i == 'C':
            r_seq = r_seq + 'G'
        elif i == 'G':
            r_seq = r_seq + 'C'
        else:
            r_seq = r_seq + i
    return r_seq


def get_End_motif(Emotif, seq1, seq2):
    if seq1.count('N') + seq1.count('n') + seq2.count('N') + seq2.count('n') != 0:
        return Emotif
    seq2 = reverse_seq(seq2)
    if seq1 in Emotif.keys():
        Emotif[seq1] = Emotif[seq1] + 1
    if seq2 in Emotif.keys():
        Emotif[seq2] = Emotif[seq2] + 1
    return Emotif


def calc_MDS(inputEndMotifFile, outputfile):
    inputfile = pd.read_table(inputEndMotifFile, header=None, names=['bases', 'frequency'])
    k_mer = math.log(len(inputfile), 4)
    frequency = inputfile['frequency'].to_numpy()
    MDS = np.sum(-frequency * np.log2(frequency) / np.log2(4 ** k_mer))
    with open(outputfile, 'a') as f:
        f.write(inputEndMotifFile + '\t' + str(MDS) + '\n')


def get_Breakpoint_motif(Bpmotif, seq1, seq2):
    # seq1 and seq2 do not include N
    if seq1.count('N') + seq1.count('n') + seq2.count('N') + seq2.count('n') != 0:
        return Bpmotif
    seq2 = reverse_seq(seq2)
    if seq1 in Bpmotif.keys():
        Bpmotif[seq1] = Bpmotif[seq1] + 1
    if seq2 in Bpmotif.keys():
        Bpmotif[seq2] = Bpmotif[seq2] + 1
    return Bpmotif


def bam_process(bamInput,
                blacklistInput,
                bedOutput,
                genome_reference,
                CHR,
                mapQuality,
                k_mer,
                fragFilter=False,
                minLen=None,
                maxLen=None):
    bedOutput_path = os.path.realpath(bedOutput)
    this_pid = os.getpid()
    tmp_split = os.path.splitext(bedOutput_path)
    tmp_bedOutput = tmp_split[0] + "-temp-" + str(this_pid) + tmp_split[1]
    bai = bamInput + '.bai'
    if not os.path.exists(bai):
        pysam.sort("-o", bamInput, bamInput)
        pysam.index(bamInput)

        message = "Index file " + bai + " do not exist!" + '  cfDNAFE use samtools to sort!'
        print(message)
    input_file = pysam.Samfile(bamInput)
    genome = pysam.Fastafile(genome_reference)
    bases = ['A', 'C', 'T', 'G']
    End_motif = {}
    End_motif_output = os.path.basename(tmp_split[0]) + '.EndMotif'
    Breakpoint_motif = {}
    Break_motif_output = os.path.basename(tmp_split[0]) + '.BreakPointMotif'
    MDS_output = os.path.basename(tmp_split[0]) + '.MDS'
    EDM_output_path = os.path.join(os.path.dirname(bedOutput_path), 'EDM/')
    BPM_output_path = os.path.join(os.path.dirname(bedOutput_path), 'BPM/')
    MDS_output_path = os.path.join(os.path.dirname(bedOutput_path), 'MDS/')
    if not os.path.exists(EDM_output_path):
        os.mkdir(EDM_output_path)
    if not os.path.exists(BPM_output_path):
        os.mkdir(BPM_output_path)
    if not os.path.exists(MDS_output_path):
        os.mkdir(MDS_output_path)
    # generated bases seq
    for i in itertools.product(*([bases] * k_mer)):
        End_motif["".join(i)] = 0
        Breakpoint_motif["".join(i)] = 0
    compress = True
    bedWrite = open(tmp_bedOutput, "w")
    print("input file:", bamInput, blacklistInput)
    print("output file:", bedOutput)
    for read1, read2 in read_pair_generator(input_file):
        if read1.mapping_quality < mapQuality or read2.mapping_quality < mapQuality or read1.reference_name not in CHR:
            continue
        read1Start = read1.reference_start
        read1End = read1.reference_end
        read2Start = read2.reference_start
        read2End = read2.reference_end
        if not read1.is_reverse:  # read1 is forward strand, read2 is reverse strand
            rstart = read1Start  # 0-based left-most site
            rend = read2End
            forward_end5 = read1.seq[0:k_mer].upper()
            forward_end3 = read2.seq[-k_mer:].upper()
        else:  # read1 is reverse strand, read2 is forward strand
            rstart = read2Start  # 0-based left-most site
            rend = read1End
            forward_end5 = read2.seq[0:k_mer].upper()
            forward_end3 = read1.seq[-k_mer:].upper()
        if (rstart < 0) or (rend < 0) or (rstart >= rend):
            continue
        if fragFilter:
            readLen = rend - rstart
            if (readLen <= minLen) or (readLen >= maxLen):
                continue
            # WPS
        gc = GCcontent(genome.fetch(read1.reference_name, rstart, rend))
        tmp_str = read1.reference_name + '\t' + str(rstart + 1) + '\t' + str(rend + 1) + '\t' + str(gc) + '\n'
        bedWrite.write(tmp_str)
        End_motif = get_End_motif(End_motif, forward_end3, forward_end3)
        pos = math.ceil(k_mer / 2)

        if k_mer % 2 == 0:
            ref_seq1 = genome.fetch(read1.reference_name, rstart - pos, rstart).upper()
            ref_seq2 = genome.fetch(read2.reference_name, rend, rend + pos).upper()
            try:
                Breakpoint_motif = get_Breakpoint_motif(Breakpoint_motif, ref_seq1 + forward_end5[0:pos],
                                                        forward_end3[-pos:] + ref_seq2)
            except:
                continue
        else:
            ref_seq1 = genome.fetch(read1.reference_name, rstart - pos + 1, rstart).upper()
            ref_seq2 = genome.fetch(read2.reference_name, rend, rend + pos - 1).upper()
            try:
                Breakpoint_motif = get_Breakpoint_motif(Breakpoint_motif, ref_seq1 + forward_end5[0:pos],
                                                        forward_end3[-pos:] + ref_seq2)
            except:
                continue

    bedWrite.close()
    print("Fragments generated and filter blackregions, waiting for sorting......")
    bedData = pybedtools.BedTool(tmp_bedOutput)
    black_reigon = pybedtools.BedTool(blacklistInput)
    bedData = bedData.subtract(black_reigon, A=True)
    bedData.sort(output=bedOutput)
    os.remove(tmp_bedOutput)
    print("Fragments sorted.")
    if compress:
        print("Waitting for compressing and indexing......")
        bedgzfile = bedOutput + ".gz"
        pysam.tabix_compress(bedOutput, bedgzfile, force=False)
        pysam.tabix_index(bedgzfile, preset="bed", zerobased=True)
        print("Indexing bedgz file finished!")

    print("motif feature generated, waiting......")

    with open(os.path.join(EDM_output_path, End_motif_output), 'w') as f:
        sum_frequency = np.sum(list(End_motif.values()))
        for k, v in End_motif.items():
            f.write(k + '\t' + str(v / sum_frequency) + '\n')
    with open(os.path.join(BPM_output_path, Break_motif_output), 'w') as f:
        sum_frequency = np.sum(list(Breakpoint_motif.values()))
        for k, v in Breakpoint_motif.items():
            f.write(k + '\t' + str(v / sum_frequency) + '\n')

    calc_MDS(os.path.join(EDM_output_path, End_motif_output), os.path.join(MDS_output_path, MDS_output))
    return True


def GCcorrect(coverage, bias):
    covl = len(coverage)
    valid = [True for i in range(covl)]
    temp_cov = []
    temp_bias = []
    for i in range(covl):
        if np.isnan(bias[i]):
            valid[i] = False
        else:
            temp_cov.append(coverage[i])
            temp_bias.append(bias[i])
    med = np.median(temp_cov)
    correct_cov = []
    i = np.arange(np.min(temp_bias), np.max(temp_bias), 0.001)

    coverage_trend = loess(temp_bias, temp_cov, span=0.75)
    coverage_trend.fit()
    coverage_model = loess(i, coverage_trend.predict(i, stderror=True).values)
    coverage_model.fit()
    coverage_pred = coverage_model.predict(temp_bias, stderror=True)
    pred = np.array(coverage_pred.values)
    coverage_corrected = temp_cov - pred + med
    i, j = 0, 0
    while i < covl:
        if valid[i]:
            if coverage_corrected[j] < 0:
                correct_cov.append(0)
            else:
                correct_cov.append(coverage_corrected[j])
            j += 1
        else:
            correct_cov.append(0)
        i += 1
    return correct_cov


def calc_FSR(bedgzInput, binInput, windows, continue_N, outputfile):
    print("input file:", bedgzInput, binInput)
    inputbed = pysam.Tabixfile(filename=bedgzInput, mode="r")
    bins = pybedtools.BedTool(binInput)
    length = len(bins)
    shorts_data, intermediates_data, longs_data, totals_data, bingc = [], [], [], [], []
    chrom = []
    print("output file:", outputfile)
    for id in range(length):
        bin = bins[id]
        try:
            chrom.append(bin.chrom)
            inputbed.fetch(bin.chrom, bin.start, bin.end)
        except ValueError:
            # no overlap read, we append value is 0
            bingc.append(np.nan)
            shorts_data.append(0)
            intermediates_data.append(0)
            longs_data.append(0)
        else:

            bin_data = []
            gc = []

            for read in inputbed.fetch(bin.chrom, bin.start, bin.end):
                bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                if 65 <= int(read.split("\t")[2]) - int(read.split("\t")[1]) <= 400:
                    gc.append(float(read.split("\t")[3]))

            count = np.bincount(bin_data, minlength=401)
            if len(gc) == 0:
                bingc.append(np.nan)
            else:
                bingc.append(np.mean(gc))

            shorts = sum(count[65: 150])
            intermediates = sum(count[151: 260])
            longs = sum(count[261: 400])
            totals = sum(count[65: 400])
            if totals == 0:
                shorts_data.append(0)
                intermediates_data.append(0)
                longs_data.append(0)
            else:
                shorts_data.append(shorts / totals)
                intermediates_data.append(intermediates / totals)
                longs_data.append(longs / totals)
            print(shorts)
    start = 0
    step = 0
    FSRfile = open(outputfile, 'w')
    FSRfile.write(
        "region" + '\t' + 'short-ratio' + '\t' + 'itermediate-ratio' + '\t' + 'long-ratio' + '\n')
    while step < length:
        num = chrom.count(chrom[step])
        continues_bin = num // continue_N
        last_bin = num % continue_N
        for i in range(continues_bin):
            tmp_array = np.zeros(3)
            bin_start = start * windows
            bin_end = (start + continue_N) * windows - 1
            combine_shorts = shorts_data[step: step + continue_N]
            combine_intermediates = intermediates_data[step: step + continue_N]
            combine_longs = longs_data[step: step + continue_N]
            tmp_array[0] = np.mean(combine_shorts)
            tmp_array[1] = np.mean(combine_intermediates)
            tmp_array[2] = np.mean(combine_longs)
            region = chrom[step] + ':' + str(bin_start) + "-" + str(bin_end)
            temp_str = region + '\t' + '\t'.join(map(str, tmp_array)) + '\n'
            FSRfile.write(temp_str)
            step += continue_N
            start += continue_N
        if last_bin != 0:
            step += last_bin
            start = 0

    FSRfile.close()
    return True


def calc_FSC(bedgzInput, binInput, windows, continue_N, outputfile):
    print("input file:", bedgzInput, binInput)
    inputbed = pysam.Tabixfile(filename=bedgzInput, mode="r")
    bins = pybedtools.BedTool(binInput)
    length = len(bins)
    shorts_data, intermediates_data, longs_data, totals_data, bingc = [], [], [], [], []
    chrom = []
    print("output file:", outputfile)
    for id in range(length):
        bin = bins[id]
        try:
            chrom.append(bin.chrom)
            inputbed.fetch(bin.chrom, bin.start, bin.end)
        except ValueError:
            # no overlap read, we append value is 0
            bingc.append(np.nan)
            shorts_data.append(0)
            intermediates_data.append(0)
            longs_data.append(0)
            totals_data.append(0)
        else:

            bin_data = []
            gc = []

            for read in inputbed.fetch(bin.chrom, bin.start, bin.end):
                bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                if 65 <= int(read.split("\t")[2]) - int(read.split("\t")[1]) <= 400:
                    gc.append(float(read.split("\t")[3]))

            count = np.bincount(bin_data, minlength=401)
            if len(gc) == 0:
                bingc.append(np.nan)
            else:
                bingc.append(np.mean(gc))
            shorts = sum(count[65: 150])
            intermediates = sum(count[151: 260])
            longs = sum(count[261: 400])
            totals = sum(count[65: 400])
            shorts_data.append(shorts)
            intermediates_data.append(intermediates)
            longs_data.append(longs)
            totals_data.append(totals)

    correct_shorts = GCcorrect(shorts_data, bingc)
    correct_intermediates = GCcorrect(intermediates_data, bingc)
    correct_longs = GCcorrect(longs_data, bingc)
    correct_totals = GCcorrect(totals_data, bingc)
    start = 0
    step = 0
    FSCfile = open(outputfile, 'w')
    short_s, intermediate_s, long_s, total_s = [], [], [], []
    region = []
    while step < length:
        num = chrom.count(chrom[step])
        continues_bin = num // continue_N
        last_bin = num % continue_N
        for i in range(continues_bin):
            bin_start = start * windows
            bin_end = (start + continue_N) * windows - 1
            combine_shorts = correct_shorts[step: step + continue_N]
            combine_intermediates = correct_intermediates[step: step + continue_N]
            combine_longs = correct_longs[step: step + continue_N]
            combine_totals = correct_totals[step: step + continue_N]
            short_s.append(np.sum(combine_shorts))
            intermediate_s.append(np.sum(combine_intermediates))
            long_s.append(np.sum(combine_longs))
            total_s.append(np.sum(combine_totals))
            region.append(chrom[step] + ':' + str(bin_start) + "-" + str(bin_end))

            step += continue_N
            start += continue_N
        if last_bin != 0:
            step += last_bin
            start = 0
    short_z = (short_s - np.mean(short_s)) / np.std(short_s)
    intermediate_z = (intermediate_s - np.mean(intermediate_s)) / np.std(intermediate_s)
    long_z = (long_s - np.mean(long_s)) / np.std(long_s)
    total_z = (total_s - np.mean(total_s)) / np.std(total_s)
    FSCfile.write(
        "region" + '\t' + 'short-fragment-zscore' + '\t' + 'itermediate-fragment-zscore' + '\t' + 'long-fragment-zscore' + '\t' + 'total-fragment-zscore' + '\n')
    for j in range(len(region)):
        temp_str = region[j] + '\t' + str(short_z[j]) + '\t' + str(intermediate_z[j]) + '\t' + str(
            long_z[j]) + '\t' + str(total_z[j]) + '\n'
        FSCfile.write(temp_str)
    FSCfile.close()
    return True


def calc_FSD(bedgzInput, armsfile, outputfile):
    print("input file:", bedgzInput, armsfile)
    inputbed = pysam.Tabixfile(filename=bedgzInput, mode="r")
    bins = pybedtools.BedTool(armsfile)
    length = len(bins)
    interval_data = []
    region = []
    print("output file:", outputfile)
    for id in range(length):
        bin = bins[id]
        region.append(bin.chrom + ':' + str(bin.start) + "-" + str(bin.end))
        try:
            inputbed.fetch(bin.chrom, bin.start, bin.end)


        except ValueError:
            # no overlap read, continue
            interval_data.append([0] * 67)
            continue
        else:
            bin_data = []
            for read in inputbed.fetch(bin.chrom, bin.start, bin.end):
                bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
            step_size = 5
            start_bin = 65
            end_bin = 400
            bin_len = int((end_bin - start_bin) / step_size)
            temp_bin = []
            count = np.bincount(bin_data, minlength=401)
            for bin_id in range(bin_len):
                temp_bin.append(count[(start_bin + step_size * bin_id):(start_bin + step_size * (bin_id + 1))])
            interval_data.append(temp_bin)

    FSDfile = open(outputfile, 'w')
    sbin = np.arange(65, 400, 5)
    head_str = 'region' + '\t'
    for s in sbin:
        head_str = head_str + str(s) + '-' + str(s+4)
    for i in range(length):
        arms = interval_data[i]
        score = np.zeros(67)

        for j in range(67):
            if np.sum(arms) != 0:
                score[j] = np.sum(arms[j]) / np.sum(arms)
        temp_str = region[i] + '\t' + '\t'.join(map(str, score)) + '\n'
        FSDfile.write(temp_str)
    FSDfile.close()
    return True


def calc_WPS(bedgzInput, tsvInput, outputfile, empty=None, protectInput=120, minSize=120, maxSize=180):
    bedgzfile = bedgzInput.strip("""\'""")
    tbx = pysam.TabixFile(bedgzfile)
    outputfile = outputfile.strip("""\'""")
    protection = protectInput // 2
    infile = open(tsvInput)
    prefix = "chr"
    validChroms = set(map(str, list(range(1, 23)) + ["X"]))  # human genome
    print("input file:", bedgzInput, tsvInput)
    for line in infile.readlines():
        (cid, chrom, start, end, strand,) = line.split()
        chrom = chrom.replace("chr", "")
        if chrom not in validChroms:
            continue
        regionStart, regionEnd = int(float(start)), int(float(end))
        if regionStart < 1:
            continue  # invalid region
        posRange = defaultdict(lambda: [0, 0])
        filteredReads = Intersecter()
        try:  # if tbx.fetch do not find any row, next row
            for row in tbx.fetch(
                    prefix + chrom, regionStart - protection, regionEnd + protection
            ):  # all fragments overlaped with this region is collected
                tmp_row = row.split()
                rstart = int(tmp_row[1])  # convert to 1-based
                rend = int(tmp_row[2])  # end included
                lseq = int(tmp_row[2]) - int(tmp_row[1])  # fragment length
                if lseq < minSize or lseq > maxSize:
                    continue
                filteredReads.add_interval(Interval(rstart, rend))  # save the fragments overlap with region
                for i in range(rstart, rend):
                    # for a single nucleotide site, compute how many reads overlaped span it (include read end point)
                    if regionStart <= i <= regionEnd:
                        posRange[i][0] += 1
                if regionStart <= rstart <= regionEnd:
                    # for a single nucleotide site, compute how many read end point located at this site
                    posRange[rstart][1] += 1
                if regionStart <= rend <= regionEnd:
                    posRange[rend][1] += 1

        except Exception:
            continue
        #
        filename = outputfile % cid
        outfile = gzip.open(filename, 'w')
        cov_sites = 0
        outLines = []
        for pos in range(regionStart, regionEnd + 1):
            rstart, rend = pos - protection, pos + protection
            gcount, bcount = 0, 0
            for read in filteredReads.find(rstart, rend):
                if (read.start > rstart) or (read.end < rend):
                    bcount += 1  # fragments located in window
                else:
                    gcount += 1  # fragments spanned window
            covCount, startCount = posRange[pos]
            cov_sites += covCount
            # chrom: chromatin, pos: position in the genome, covCount:how many reads span this site, startCount: how many reads end point located
            # in this site, gcount-bcount: WPS
            outLines.append("%s\t%d\t%d\t%d\t%d\n" % (chrom, pos, covCount, startCount, gcount - bcount))
        if strand == "-":
            outLines = outLines[::-1]  # - strand!!!
        for line in outLines:
            outfile.write(line.encode())  # write in binary
        outfile.close()
        if cov_sites == 0 and not empty:  # remove empty files
            os.remove(filename)
    return True


def calc_OCF(bedgzInput, ocrInput, outputdir):
    print("input file:", bedgzInput, ocrInput)
    tbx = pysam.TabixFile(bedgzInput)
    regions = pd.read_csv(ocrInput, sep="\t", header=None, names=["chr", "start", "end", "description"])

    leftPOS = defaultdict(partial(defaultdict, int))
    rightPOS = defaultdict(partial(defaultdict, int))
    total = defaultdict(lambda: [0, 0])
    for idx, region in regions.iterrows():
        region_Chr, region_Start, region_End, region_Label = (
            region["chr"],
            region["start"],
            region["end"],
            region["description"]
        )

        if region_Start < 1:
            message = "Start of the region must > 0!"

        # fetch read in the region
        try:
            fetched_reads = tbx.fetch(region_Chr, region_Start, region_End)
        except ValueError:
            continue
        for row in fetched_reads:
            tmp_row = row.split()
            rstart = int(tmp_row[1])  # convert to 1-based
            rend = int(tmp_row[2])  # end included

            if rstart >= region_Start:
                s = rstart - region_Start
                leftPOS[region_Label][s] += 1
                total[region_Label][0] += 1
            if rend <= region_End:
                e = rend - region_Start + 1
                rightPOS[region_Label][e] += 1
                total[region_Label][1] += 1
    Labels = []
    ocf = []
    outputfile = os.path.join(outputdir, 'all.ocf.csv')
    print("output file: ", outputfile)
    for label in total.keys():
        output = os.path.join(outputdir, label) + '.sync.end'
        output_write = open(output, 'w')
        Labels.append(label)
        le = leftPOS[label]
        re = rightPOS[label]
        ts = total[label][0] / 10000
        te = total[label][1] / 10000
        num = 2000
        for k in range(num):
            l = le[k]
            r = re[k]
            output_write.write(
                str(k - 1000) + '\t' + str(l) + '\t' + str(l / ts) + '\t' + str(r) + '\t' + str(r / te) + '\n')
        output_write.close()
        with open(output, 'r') as o:
            peak = 60
            bin = 10
            trueends = 0  # neg (-peak-bin to -peak+bin) + pos (peak-bin to peak+bin)
            background = 0  # neg (peak-bin to peak+bin) + pos (-peak-bin to -peak+bin)
            for line in o.readlines():
                (loc, left, Left, right, Right,) = line.split()
                if -peak - bin <= int(loc) <= -peak + bin:
                    trueends += float(Right)
                    background += float(Left)
                elif peak - bin <= int(loc) <= peak + bin:
                    trueends += float(Left)
                    background += float(Right)
                else:
                    continue
            ocf.append(trueends - background)
        ocf_df = pd.DataFrame({"tissue": Labels, "OCF": ocf})
        ocf_df.to_csv(outputfile, sep="\t", index=None)


def calc_UXM(bamInput, markInput, outputFile, mapQuality, minCpG, methyThreshold, unmethyThreshold, type):
    bai = bamInput + '.bai'
    if not os.path.exists(bai):
        pysam.sort("-o", bamInput, bamInput)
        pysam.index(bamInput)
        message = "Index file " + bai + " do not exist!" + '  cfDNAFE use samtools to sort!'
        print(message)
    input_file = pysam.Samfile(bamInput)
    marks = pybedtools.BedTool(markInput)
    print("input file:", bamInput, markInput)
    print("output file:", outputFile)
    res = []
    for i in range(len(marks)):
        mark = marks[i]
        region = mark.chrom + ':' + str(mark.start) + "-" + str(mark.end)
        try:
            input_file.fetch(mark.chrom, mark.start, mark.end)
        except ValueError:
            # no overlap read, we append value is 0
            res.append(region + "\t" + '0' + '\t' + '0' + '\t' + '0')
        else:
            Ufragment = 0  # mostly unmethylated fragment
            Xfragment = 0  # mixed fragment
            Mfragment = 0  # mostly methylated fragment
            if type == 'PE':
                for read1, read2 in read_pair_generator(input_file, region):
                    if read1.mapping_quality < mapQuality or read2.mapping_quality < mapQuality:
                        continue
                    read1Start = read1.reference_start
                    read1End = read1.reference_end
                    read2Start = read2.reference_start
                    read2End = read2.reference_end
                    if not read1.is_reverse:  # read1 is forward strand, read2 is reverse strand
                        if (read2Start < read1End):
                            overlap = read1End - read2Start
                            m1 = read1.get_tag("XM")
                            m2 = read2.get_tag("XM")
                            num_methylated = m1.count("Z") + m2[overlap:].count("Z")
                            num_unmethylated = m1.count("z") + m2[overlap:].count("z")
                        else:
                            m1 = read1.get_tag("XM")
                            m2 = read2.get_tag("XM")
                            num_methylated = m1.count("Z") + m2.count("Z")
                            num_unmethylated = m1.count("z") + m2.count("z")

                    else:  # read1 is reverse strand, read2 is forward strand
                        if (read1Start < read2End):
                            overlap = read2End - read1Start
                            m1 = read2.get_tag("XM")
                            m2 = read1.get_tag("XM")
                            num_methylated = m2.count("Z") + m1[overlap:].count("Z")
                            num_unmethylated = m2.count("z") + m1[overlap:].count("z")
                        else:
                            m1 = read1.get_tag("XM")
                            m2 = read2.get_tag("XM")
                            num_methylated = m1.count("Z") + m2.count("Z")
                            num_unmethylated = m1.count("z") + m2.count("z")

                    if num_methylated + num_unmethylated < minCpG:
                        continue

                    if num_methylated / (num_methylated + num_unmethylated) >= methyThreshold:
                        Mfragment = Mfragment + 1
                    elif num_methylated / (num_methylated + num_unmethylated) <= unmethyThreshold:
                        Ufragment = Ufragment + 1
                    else:
                        Xfragment = Xfragment + 1

            elif type == 'SE':
                for read in input_file.fetch(mark.chrom, mark.start, mark.end):
                    if read.mapping_quality < mapQuality:
                        continue
                    m = read.get_tag("XM")
                    num_methylated = m.count("Z")
                    num_unmethylated = m.count("z")
                    if num_methylated + num_unmethylated < minCpG:
                        continue

                    if num_methylated / (num_methylated + num_unmethylated) >= methyThreshold:
                        Mfragment = Mfragment + 1
                    elif num_methylated / (num_methylated + num_unmethylated) <= unmethyThreshold:
                        Ufragment = Ufragment + 1
                    else:
                        Xfragment = Xfragment + 1
            else:
                raise commonError("Parameter type must be SE or PE")
            total = Mfragment + Ufragment + Xfragment
            if total == 0:
                res.append(region + '\t' + '0' + '\t' + '0' + '\t' + '0')
            else:
                tmp_array = np.zeros(3)

                tmp_array[0] = Ufragment / total
                tmp_array[1] = Xfragment / total
                tmp_array[2] = Mfragment / total

                res.append(region + '\t' + '\t'.join(map(str, tmp_array)))

    with open(outputFile, 'w') as f:
        f.write('region' + '\t' + 'U' + '\t' + 'X' + '\t' + 'M' + '\n')
        for i in res:
            f.write(i + '\n')
