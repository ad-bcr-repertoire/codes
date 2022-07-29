import csv
import os
import re
from Common.SeqUtil_test import Helper
from Common.ProcessUtil import run_clustal_omega
import glob
from shutil import copyfile
import time
from Common.FileUtil import read_file, write_file_uncond
from Common.Enum import ChainType

IGHC_REFERENCE_FILE = '/home/Human/reference/IGHC_CH1.csv'

forwardPrimerList = ['GGCCTCAGTGAAGGTCTCCTGCAAG', 'GTCTGGTCCTACGCTGGTGAACCC', 'CTGGGGGGTCCCTGAGACTCTCCTG',
                     'CTTCGGAGACCCTGTCCCTCACCTG', 'CGGGGAGTCTCTGAAGATCTCCTGT', 'TCGCAGACCCTCTCACTCACCTGTG']
reversePrimerList = ['GAAGGAAGTCCTGTGCGAG', 'GGGAAGTAGTCCTTGACCA', 'GGGGAAGAAGCCCTGGAC',
                    'TGGGTGGTACCCAGTTATCAA', 'AAGTAGCCCGTGGCCAGG']

def read_fasta(inputFile):
    resultSeqList = []
    with open(inputFile, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                if line[0] != '>':
                    print('File format error!')
                    raise RuntimeError
                tempSeq = {'name':line[1:-1], 'sequence':''}
            elif line[0] == '>':
                resultSeqList.append(tempSeq)
                tempSeq = {'name':line[1:-1], 'sequence':''}
            else:
                tempSeq['sequence'] += line[:-1]
        resultSeqList.append(tempSeq)
    return resultSeqList

def primer_recognition_and_UMI_extraction(input_file_name, output_dir, forward_primer_list, reverse_primer_list):
    sHelper = Helper()

    reversePrimerList = [sHelper.to_reverse_complement(reversePrimer) for reversePrimer in reverse_primer_list]

    reversePrimerMatchList = []
    for reversePrimer in reversePrimerList:
        reversePrimerMatchList.append(reversePrimer)
    for reversePrimer in reversePrimerList:
        for i in range(len(reversePrimer)):
            reversePrimerMatchList.append(reversePrimer[:i] + reversePrimer[i + 1:])
            reversePrimerMatchList.append(reversePrimer[:i] + '[ATGC]' + reversePrimer[i + 1:])

    forwardPrimerMatchList = []
    for forwardPrimer in forward_primer_list:
        forwardPrimerMatchList.append(forwardPrimer)
    for forwardPrimer in forward_primer_list:
        for i in range(len(forwardPrimer)):
            forwardPrimerMatchList.append(forwardPrimer[:i] + forwardPrimer[i + 1:])
            forwardPrimerMatchList.append(forwardPrimer[:i] + '[ATGC]' + forwardPrimer[i + 1:])

    reverseTrimNum = 40
    forwardTrimNum = 30

    newSeqList = []

    UMISet = {}

    with open(os.path.join(input_file_name), 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

        if len(header) != 2:
            print('Input file column error')
            raise RuntimeError

        for row in reader:
            sequence = row[0]
            readcount = int(row[1])

            for reversePrimerMatch in reversePrimerMatchList:
                reverseMatchResult = re.search(reversePrimerMatch, sequence[-reverseTrimNum:])

                if reverseMatchResult != None:
                    UMI = sequence[-(reverseTrimNum - reverseMatchResult.span()[1]):]

                    for forwardPrimerMatch in forwardPrimerMatchList:
                        forwardMatchResult = re.search(forwardPrimerMatch, sequence[:forwardTrimNum])

                        if forwardMatchResult != None:
                            newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[0])]
                            newSeqList.append({'UMI':UMI, 'full_NT':newSequence, 'readcount':readcount})

                            if len(UMI) == 14 or len(UMI) == 15 or len(UMI) == 13:
                                try:
                                    UMISet[UMI].append({'UMI':UMI, 'full_NT':newSequence, 'readcount':readcount})
                                except KeyError:
                                    UMISet[UMI] = [{'UMI':UMI, 'full_NT':newSequence, 'readcount':readcount}]
                            break
                    break

    rcUMISet = {}

    for UMI in UMISet:
        readcountSum = 0
        for each in UMISet[UMI]:
            readcountSum += int(each['readcount'])
        try:
            rcUMISet[str(readcountSum)].append(UMI)
        except KeyError:
            rcUMISet[str(readcountSum)] = [UMI]

    rcList = [int(x) for x in list(set(rcUMISet.keys()))]
    rcList.sort(reverse=True)

    mergedUMISet = {}

    sortedUMISet = {}

    for rc in rcList:
        for rcUMI in rcUMISet[str(rc)]:
            sortedUMISet[rcUMI] = UMISet[rcUMI]

    while (len(sortedUMISet) > 0):
        sortedUMIList = list(sortedUMISet.keys())
        queryUMI = sortedUMIList.pop(0)
        del sortedUMISet[queryUMI]
        UMIMatchList = []

        woNTSet = {'A':['T','G','C'], 'T':['A','G','C'], 'G':['A','T','C'], 'C':['A','T','G']}

        for i in range(len(queryUMI)):
            UMIMatchList.append(queryUMI[:i] + queryUMI[i + 1:])
            for woNT in woNTSet[queryUMI[i]]:
                UMIMatchList.append(queryUMI[:i] + woNT + queryUMI[i + 1:])

        UMIMatchList = list(set(UMIMatchList))

        toBeMergedUMIList = []

        for UMIMatch in UMIMatchList:
            try:
                sortedUMISet[UMIMatch]
                toBeMergedUMIList.append(UMIMatch)
            except KeyError:
                pass

        if len(toBeMergedUMIList) == 0:
            mergedUMISet[queryUMI] = UMISet[queryUMI]
        else:
            mergedUMISet[queryUMI] = UMISet[queryUMI]
            for toBeMergedUMI in toBeMergedUMIList:
                mergedUMISet[queryUMI] += UMISet[toBeMergedUMI]
                del sortedUMISet[toBeMergedUMI]

    for mergedUMI in mergedUMISet:
        mergedSeqNum = len(mergedUMISet[mergedUMI])
        mergedRcSum = 0
        for mergedUMISeq in mergedUMISet[mergedUMI]:
            mergedRcSum += mergedUMISeq['readcount']

        with open(os.path.join(output_dir, mergedUMI + '_' + str(mergedSeqNum) + '_' + str(mergedRcSum) + '.fasta'), 'w') as f:
            for mergedUMISeq in mergedUMISet[mergedUMI]:
                f.write('>' + str(mergedUMISeq['readcount']) + '\n')
                f.write(mergedUMISeq['full_NT'] + '\n')


def sub_cluster_UMI_data(input_dir, output_dir, dist_thresh):
    for inputFile in glob.glob(os.path.join(input_dir, '*.fasta')):
        UMI = os.path.basename(inputFile).split('.')[0].split('_')[0]
        readcountSet = {}

        fastaList = read_fasta(inputFile)
        for fasta in fastaList:
            tempSeq = {'full_NT':fasta['sequence'], 'readcount':int(fasta['name'])}
            try:
                readcountSet[fasta['name']].append(tempSeq)
            except KeyError:
                readcountSet[fasta['name']] = [tempSeq]

        readcountList = []
        for readcount in readcountSet:
            readcountList.append(int(readcount))

        readcountList = list(set(readcountList))
        readcountList.sort(reverse=True)

        seqList = []
        for readcount in readcountList:
            seqList += readcountSet[str(readcount)]

        subClusterNum = 1
        while(len(seqList) > 0):
            targetSeq = seqList.pop(0)

            resultSeqList = [targetSeq]
            readcountSum = targetSeq['readcount']

            popIndexList = []

            for i, seq in enumerate(seqList):
                distFromTarget = Helper.levenshtein_distance(targetSeq['full_NT'], seq['full_NT'])
                if distFromTarget <= dist_thresh:
                    resultSeqList.append(seq)
                    readcountSum += seq['readcount']
                    popIndexList.append(i)

            popIndexList.sort(reverse=True)

            for popIndex in popIndexList:
                seqList.pop(popIndex)

            if readcountSum > 1:
                with open(os.path.join(output_dir, UMI + '-' + str(subClusterNum) + '_' + str(len(resultSeqList)) + '_' + str(readcountSum) + '.fasta'), 'w') as f:
                    for resultSeq in resultSeqList:
                        f.write('>' + str(resultSeq['readcount']) + '\n')
                        f.write(resultSeq['full_NT'] + '\n')
                subClusterNum += 1

def msa_by_clustal_omega(input_dir, output_dir):
    for inputFile in glob.glob(os.path.join(input_dir, '*.fasta')):
        inputFileName = os.path.basename(inputFile)
        inputFileNameSplit = inputFileName.split('.')[0].split('_')
        seqNum = int(inputFileNameSplit[1])

        if seqNum > 1:
            run_clustal_omega(input=inputFile, output=os.path.join(output_dir, inputFileName))
        else:
            copyfile(inputFile, os.path.join(CLUSTAL_OMEGA_DIR, inputFileName))

def extract_consensus_sequence(input_dir, vote_thresh, output_file):
    resultSeqSet = {}
    for inputFile in glob.glob(os.path.join(input_dir, '*.fasta')):
        inputFileName = os.path.basename(inputFile)
        inputFileNameSplit = inputFileName.split('.')[0].split('_')
        seqNum = int(inputFileNameSplit[1])
        rcSum = int(inputFileNameSplit[2])

        if rcSum == 1:
            continue

        if seqNum == 2 and rcSum == 2:
            continue

        seqList = read_fasta(inputFile)

        seqScoreList = []

        for i, seq in enumerate(seqList):
            if i == 0:
                for eachSeq in seq['sequence']:
                    tempSet = {'A':0, 'G':0, 'T':0, 'C':0, '-':0}
                    tempSet[eachSeq] += int(seq['name'])
                    seqScoreList.append(tempSet)
            else:
                for j, eachSeq in enumerate(seq['sequence']):
                    seqScoreList[j][eachSeq] += int(seq['name'])

        resultSeq = ''
        for seqScore in seqScoreList:
            maxVote = max(seqScore.values())
            if float(maxVote)/rcSum > vote_thresh:
                maxVoteNT = list(seqScore.keys())[list(seqScore.values()).index(maxVote)]
                resultSeq += maxVoteNT
            else:
                resultSeq += 'X'

        try:
            resultSeq.index('X')
        except ValueError:
            try:
                resultSeqSet[resultSeq.replace('-','')] += 1
            except KeyError:
                resultSeqSet[resultSeq.replace('-','')] = 1

    with open(output_file, 'w') as f:
        f.write('full_NT,readcount\n')
        for resultSeq in resultSeqSet:
            f.write(resultSeq + ',' + str(resultSeqSet[resultSeq]) +'\n')

def c_gene_annotation(input_file, output_file, blast_dir='', num_threads=10):
    seqHelper = Helper()
    chRefHeder, chRefList = read_file(in_file=IGHC_REFERENCE_FILE)
    targetSetHeader, targetSeqList = read_file(input_file)

    if blast_dir == '':
        blast_dir = os.path.dirname(output_file) + 'temp_cblast'

    if os.path.exists(blast_dir) != True:
        os.makedirs(blast_dir)

    blastResultList = seqHelper.blast_run(query_list=[x['full_NT'] for x in targetSeqList],
                                          db_list=[x['full_NT'] for x in chRefList], m_type='NT', num_alignments=1,
                                          num_threads=num_threads, blast_dir=blast_dir)

    newSeqSet = {}

    for i, blastResult in enumerate(blastResultList):
        if blastResult['num_alignments'] == 0:
            pass
        else:
            alignments = blastResult['alignments']
            alignment = alignments[0]

            if int(alignment['whole']) < 74 or int(alignment['whole']) > 94:
                continue

            if int(alignment['num_mismatches']) > 5:
                continue

            if int(alignment['sbjct_start']) > 10:
                continue

            if int(alignment['query_start']) > 400:
                continue

            sbjctName = chRefList[int(alignment['sbjct_index'])]['name']
            querySeq = targetSeqList[i]['full_NT']
            queryRC = targetSeqList[i]['readcount']

            newSeq = querySeq[:alignment['query_start'] - alignment['sbjct_start']]
            newSeqIsotype = sbjctName.split('*')[0][3:]
            newSeqRC = queryRC

            try:
                newSeqSet[newSeq + '|' + newSeqIsotype] += newSeqRC
            except KeyError:
                newSeqSet[newSeq + '|' + newSeqIsotype]  = newSeqRC

    newHeader = ['full_NT', 'isotype', 'readcount']
    newSeqList = []
    for newSeq in newSeqSet:
        seq, isotype = newSeq.split('|')
        newSeqList.append({'full_NT':seq, 'isotype':isotype, 'readcount':newSeqSet[newSeq]})
    write_file_uncond(out_file=output_file, header=newHeader, data=newSeqList)


def annotate_regions(input_file, output_file, annotation_list, blast_dir = ''):
    if blast_dir == '':
        blast_dir = os.path.join(os.path.dirname(output_file), 'temp_igblast')

    if os.path.exists(blast_dir) == False:
        os.makedirs(blast_dir)

    seqHelper = Helper()

    targetSeqList = []

    seqHeader, seqList = read_file(input_file)

    fullNTIndex = seqHeader.index('full_NT')

    inputFileName = os.path.splitext(os.path.basename(input_file))[0]

    for seq in seqList:
        targetSeqList.append(seq['full_NT'])

    resultList = seqHelper.annotate_regions(sequence_list=targetSeqList,
                                            chain_type=ChainType.HUMAN_HEAVY,
                                            file_blast=os.path.join(blast_dir, inputFileName + '_igblasted.txt'),
                                            annotation_list=annotation_list)

    for seq, result in zip(seqList, resultList):
        result.update(seq)

    resultHeader = seqHeader[:fullNTIndex + 1] + annotation_list + seqHeader[fullNTIndex + 1:]

    write_file_uncond(out_file=output_file, header=resultHeader, data=resultList)

def trim_and_merge_sequence(input_file, ouput_file, annotate_full_AA = True):
    seqHelper = Helper()
    seqHeader, seqList = read_file(input_file)
    newSeqSet = {}
    fullNTIndex = seqHeader.index('full_NT')
    newSeqHeader = seqHeader[:fullNTIndex + 1] + ['full_AA'] + seqHeader[fullNTIndex + 1:]

    for seq in seqList:
        seqIsotype = seq['isotype']
        if seq['FR1_to'] == 'N/A':
            continue
        newSeq = seq['full_NT'][seq['FR1_to'] - int(seq['FR1_to'] / 3) * 3:]
        seq['full_NT'] = newSeq
        try:
            newSeqSet[newSeq + '|' + seqIsotype]['readcount'] += seq['readcount']
        except KeyError:
            newSeqSet[newSeq + '|' + seqIsotype] = seq

    newSeqList = []
    for newSeq in newSeqSet:
        targetSeq = newSeqSet[newSeq]
        if annotate_full_AA:
            targetSeq['full_AA'] = seqHelper.translate(targetSeq['full_NT'])
        newSeqList.append(targetSeq)

    write_file_uncond(out_file=ouput_file, data=newSeqList, header=newSeqHeader)

def filter_non_functional(input_file, output_file, target_isotype_list):
    seqHeader, seqList = read_file(in_file=input_file)

    total = 0
    passed = 0

    totalRc = 0
    passedRc = 0

    newSeqList = []

    for seq in seqList:
        total += 1
        totalRc += int(seq['readcount'])

        # isotype filtering
        try:
            target_isotype_list.index(seq['isotype'])
        except ValueError:
            continue

        # gene annotation filtering
        if seq['V_gene'] == 'N/A':
            continue

        if seq['J_gene'] == 'N/A':
            continue

        # no CDRs extraction filtering
        if seq['CDR1_AA'] == 'N/A' or seq['CDR2_AA'] == 'N/A' or seq['CDR3_AA'] == 'N/A':
            continue

        # full AA filtering
        try:
            seq['full_AA'].index('X')
            continue
        except ValueError:
            pass

        try:
            seq['full_AA'].index('*')
            continue
        except ValueError:
            pass

        # CDR1 AA filtering
        try:
            seq['CDR1_AA'].index('*')
            continue
        except ValueError:
            pass

        try:
            seq['CDR1_AA'].index('X')
            continue
        except ValueError:
            pass

        # CDR2_AA filtering
        try:
            seq['CDR2_AA'].index('*')
            continue
        except ValueError:
            pass

        try:
            seq['CDR2_AA'].index('X')
            continue
        except ValueError:
            pass

        # CDR3_AA filtering
        try:
            seq['CDR3_AA'].index('*')
            continue
        except ValueError:
            pass

        try:
            seq['CDR3_AA'].index('X')
            continue
        except ValueError:
            pass

        newSeqList.append(seq)

        passed += 1
        passedRc += int(seq['readcount'])

    write_file_uncond(out_file=output_file, data=newSeqList, header=seqHeader)

    print(str(passedRc) + ' (' + str(round(float(passedRc) / totalRc * 100, 2)) + '%) reads passed')
    print(str(passed) + ' (' + str(round(float(passed)/total*100,2)) + '%) unique sequences passed')

SAMPLE_NAME = 'SU0337'

WORK_DIR = '/home/Human/data/AD/' + SAMPLE_NAME

DB_EXTRACTION_DIR = os.path.join(WORK_DIR, '1_db_extraction')
ERROR_CORRECTION_DIR = os.path.join(WORK_DIR, '2_error_correction')
ANNOTATION_DIR = os.path.join(WORK_DIR, '3_annotation')
FUNCTIONALITY_DIR = os.path.join(WORK_DIR, '4_functionality')

UMI_CLUSTER_DIR = os.path.join(DB_EXTRACTION_DIR, '1_UMI_cluster')
UMI_SUB_CLUSTER_DIR = os.path.join(DB_EXTRACTION_DIR, '2_UMI_sub_cluster')
CLUSTAL_OMEGA_DIR = os.path.join(DB_EXTRACTION_DIR, '3_clustal_omega')

if not os.path.exists(DB_EXTRACTION_DIR):
    os.makedirs(DB_EXTRACTION_DIR)

if not os.path.exists(ERROR_CORRECTION_DIR):
    os.makedirs(ERROR_CORRECTION_DIR)

if not os.path.exists(ANNOTATION_DIR):
    os.makedirs(ANNOTATION_DIR)

if not os.path.exists(FUNCTIONALITY_DIR):
    os.makedirs(FUNCTIONALITY_DIR)


if not os.path.exists(UMI_CLUSTER_DIR):
    os.makedirs(UMI_CLUSTER_DIR)

if not os.path.exists(UMI_SUB_CLUSTER_DIR):
    os.makedirs(UMI_SUB_CLUSTER_DIR)

if not os.path.exists(CLUSTAL_OMEGA_DIR):
    os.makedirs(CLUSTAL_OMEGA_DIR)


READ_FILE_NAME = SAMPLE_NAME + '.csv'

inputFileName = os.path.join(WORK_DIR, READ_FILE_NAME)

# 1. recognize primer sequence and extract UMIs (UMI merged by 1 mistmatch)
if True:
    sTime = time.time()
    primer_recognition_and_UMI_extraction(input_file_name=inputFileName, output_dir=UMI_CLUSTER_DIR, forward_primer_list=forwardPrimerList, reverse_primer_list=reversePrimerList)
    eTime = time.time()
    print(str(eTime-sTime) + 'seconds elapsed for [1] primer_recognition_and_UMI_extraction')

# 2. make sub-cluster
if True:
    sTime = time.time()
    sub_cluster_UMI_data(input_dir=UMI_CLUSTER_DIR, output_dir=UMI_SUB_CLUSTER_DIR, dist_thresh=5)
    eTime = time.time()
    print(str(eTime-sTime) + 'seconds elapsed for [2] sub_cluster_UMI_data')

# 3. run clustal omega
if True:
    sTime = time.time()
    msa_by_clustal_omega(input_dir=UMI_SUB_CLUSTER_DIR, output_dir=CLUSTAL_OMEGA_DIR)
    eTime = time.time()
    print(str(eTime-sTime) + 'seconds elapsed for [3] msa_by_clustal_omega')

# 4. extract consensus sequences by using clustal omerga results
if True:
    sTime = time.time()
    extract_consensus_sequence(input_dir=CLUSTAL_OMEGA_DIR, output_file=os.path.join(ERROR_CORRECTION_DIR, SAMPLE_NAME + '_d1_e1.csv'), vote_thresh=0.6)
    eTime = time.time()
    print(str(eTime-sTime) + 'seconds elapsed for [4] extract_consensus_sequence')

# 5. annotate C gene
if True:
    sTime = time.time()
    c_gene_annotation(input_file=os.path.join(ERROR_CORRECTION_DIR, SAMPLE_NAME + '_d1_e1.csv'), output_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'))
    eTime = time.time()
    print(str(eTime - sTime) + 'seconds elapsed for [5] c_gene_annotation')

# 6. trim non-functional region nucleotide & annotate additional info.
if True:
    annotationList = ['CDR1_NT', 'CDR2_NT', 'CDR3_NT', 'CDR1_AA', 'CDR2_AA', 'CDR3_AA',
                      'V_gene', 'D_gene', 'J_gene', 'FR1_to']

    annotate_regions(input_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'), output_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'),
                     annotation_list=annotationList)

    trim_and_merge_sequence(input_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'), ouput_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'))

# filter non-functional reads
if True:
    filter_non_functional(input_file=os.path.join(ANNOTATION_DIR, SAMPLE_NAME + '_d1_e1_a1.csv'), output_file=os.path.join(FUNCTIONALITY_DIR, SAMPLE_NAME + '_d1_e1_a1_f2.csv'),
                          target_isotype_list=['M', 'D', 'E', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2'])




