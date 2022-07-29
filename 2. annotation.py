import logging
import re
import os
from collections import defaultdict
from datetime import datetime
from multiprocessing import Process, Manager
import glob
import numpy as np
from Common.FileUtil import read_file, write_file, write_file_uncond
import random

from Common.FileUtil import read_file, write_file, write_file_uncond

from Common.SeqUtil_test import Helper
from typing import List, Dict, Tuple
import csv

from Common.Enum import ChainType

def correct_error(file_in: str, file_out: str, evalue: float, col_sample: int, col_num: int, col_seq: str = 'full_NT'):
   header, data = read_file(in_file=file_in)
   samples = header[(col_sample - 1):(col_sample + col_num - 1)]
   # print(samples)
   # Manipulating singleton-removed data
   logging.info('(2a) Singleton removal started : %s' % file_in)
   sr_data = []
   for d in data:
       readcount = 0
       for sample in samples:
           readcount += d[sample]
       if readcount > 1:
           sr_data.append(d)
   logging.info('(2a) Singleton removal finished')
   logging.info('(2b) Correct error started : %s' % file_in)
   sub_map = {
       'A': 'GCTN',
       'G': 'ACTN',
       'C': 'AGTN',
       'T': 'AGCN'
   }
   # generate dict
   map_seq_data = {}
   for d in sr_data:
       map_seq_data[d[col_seq]] = d
   remove_clones = set()
   logging.info('---- original %d clones' % len(sr_data))
   for sample in header[(col_sample - 1): (col_sample + col_num - 1)]:
       logging.info('---- correct_error sample: %s' % sample)
       col_reads = sample
       num_remove = 0
       for d in sr_data:
           seq = d[col_seq]
           reads1 = d[col_reads]
           if reads1 * evalue < 2:
               pass
           else:
               # substitution
               for i in range(len(seq)):
                   for mut_base in sub_map[seq[i]]:
                       mut_seq = seq[:i] + mut_base + seq[i + 1:]
                       if mut_seq in map_seq_data:
                           reads2 = map_seq_data[mut_seq][col_reads]
                           if reads1 * evalue > reads2:
                               # remove_clones.add(mut_seq)
                               map_seq_data[mut_seq][col_reads] = 0
                               num_remove += 1
       logging.info('---- read count correct %d clones' % num_remove)
   # remove [0, 0, 0, 0 ... ]
   clone_filter_out = []
   clone_filter_out.append([0 for number in range(col_num)])
   for j in range(col_num):
       clone_filter_out.append([0 for number in range(j)] + [1] + [0 for number in range(col_num - j - 1)])
   for seq in map_seq_data:
       reads_list = []
       for i in range((col_sample - 1), (col_sample + col_num - 1)):
           reads_list += [map_seq_data[seq][header[i]]]
       # if reads_list == clone_filter_out:
       if reads_list in clone_filter_out:
           remove_clones.add(seq)
   logging.info('---- remove zeros %d clones' % len(remove_clones))
   data_out = [d for d in sr_data if d[col_seq] not in remove_clones]
   write_file(file_out, header, data_out)
   logging.info('(2b) Correct error finished')

def set_id(file_in: str, file_out: str, col_id: str, col_seq: str, id_prefix: str):
    logging.info('(4a) Setting ID started')
    header, data = read_file(file_in)
    #data.sort(key=lambda x: x['X-PI'], reverse=True)
    map_length_rank = defaultdict(int)
    for d in data:
        rank = map_length_rank[len(d[col_seq])]
        d[col_id] = id_prefix + '%d-%d' % (len(d[col_seq]), rank + 1)
        map_length_rank[len(d[col_seq])] = rank + 1
    if col_id not in header:
        header = [col_id] + header

    write_file(file_out, header, data)

    logging.info('(4a) Setting ID finished')

class EdgeExtracterLD(Process):
    def __init__(self, output_dir: str, distance_threshold: int, col_id: str, col_seq: str, map_length_data: str,
                 in_queue):
        Process.__init__(self)
        self.in_queue = in_queue

        self.output_dir = output_dir
        self.dist_cutoff = distance_threshold
        self.col_id = col_id
        self.col_seq = col_seq
        self.map_length_data = map_length_data

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            item = self.in_queue.get()
            if item is None:
                self.in_queue.task_done()
                break

            self.do_work(item)
            self.in_queue.task_done()

    def do_work(self, length: int):
        cnt_compare, cnt_edge = 0, 0
        temp_start = datetime.now()

        data_l = self.map_length_data[length]
        line_edges = []
        # brute-force
        for i, d_from in enumerate(data_l):
            for d_to in data_l[i + 1:]:
                dist = Helper.levenshtein_distance(d_from[self.col_seq], d_to[self.col_seq])
                cnt_compare += 1
                if dist <= self.dist_cutoff:
                    if d_from[self.col_id] < d_to[self.col_id]:
                        line_edges.append('%s,%s,%d\n' % (d_from[self.col_id], d_to[self.col_id], dist))
                    else:
                        line_edges.append('%s,%s,%d\n' % (d_to[self.col_id], d_from[self.col_id], dist))
                    cnt_edge += 1

        for length_to in range(length + 3, length + self.dist_cutoff + 1, 3):
            if length_to in self.map_length_data:
                for d_to in self.map_length_data[length_to]:
                    dist = Helper.levenshtein_distance(d_from[self.col_seq], d_to[self.col_seq])
                    cnt_compare += 1
                    if dist <= self.dist_cutoff:
                        if d_from[self.col_id] < d_to[self.col_id]:
                            line_edges.append('%s,%s,%d\n' % (d_from[self.col_id], d_to[self.col_id], dist))
                        else:
                            line_edges.append('%s,%s,%d\n' % (d_to[self.col_id], d_from[self.col_id], dist))
                        cnt_edge += 1

        print("%s do_work ---------- takes %s, length %d, %d seqs, %d comps, %d edges (Brute-force)" % (
            datetime.now(), datetime.now() - temp_start, length, len(data_l), cnt_compare, cnt_edge))

        if len(line_edges) > 0:
            with open(os.path.join(self.output_dir, 'edges.length-%d.csv' % length), 'w') as f:
                f.writelines(line_edges)

def extract_edge_all_levenstein_under_thres(file_in_nodes: str, file_out_edges: str, col_id: str, col_seq: str,
                                            threshold: int = 30, num_process: int = 4):
    logging.info('Started extract_edge_levenstein_all_under_thres: %s' % os.path.basename(file_in_nodes))
    time_start = datetime.now()

    temp_start = datetime.now()
    header, data = read_file(file_in_nodes)
    logging.info("---------- takes %s : read_file(file_in_nodes)" % (datetime.now() - temp_start))

    temp_start = datetime.now()
    map_length_data = {}
    for d in data:
        if len(d[col_seq]) not in map_length_data:
            map_length_data[len(d[col_seq])] = []
        map_length_data[len(d[col_seq])].append(d)
    logging.info("---------- takes %s : make map_length_data" % (datetime.now() - temp_start))

    dir_edges = os.path.join(os.path.dirname(file_out_edges), os.path.basename(file_out_edges).replace('.csv', ''))
    if not os.path.exists(dir_edges):
        os.makedirs(dir_edges)

    # Make another worker thread to send query
    manager = Manager()
    in_queue = manager.Queue()
    worker_list = []
    logging.info("Spawning %d EdgeExtracterLD processes.", num_process)
    for i in range(num_process):
        worker = EdgeExtracterLD(dir_edges, threshold, col_id, col_seq, map_length_data, in_queue)
        worker.daemon = True
        worker.start()

        worker_list.append(worker)

    jobs = list(map_length_data.items())
    jobs.sort(key=lambda kv: len(kv[1]), reverse=True)

    for length, data_l in jobs:
        in_queue.put_nowait(length)

    # enqueue stopper
    for i in range(len(worker_list)):
        in_queue.put_nowait(None)

    for worker in worker_list:
        worker.join()

    logging.info('Start merging')

    edges_files = glob.glob(os.path.join(dir_edges, '*.csv'))
    with open(file_out_edges, 'w') as fw:
        fw.write('from,to,dist\n')
        for fr in edges_files:
            with open(fr, 'r') as hr:
                for line in hr:
                    fw.write(line)

    logging.info('Finished extract_edge_levenstein_all_under_thres: %s' % os.path.basename(file_in_nodes))
    logging.info('------------- takes total %s' % (datetime.now() - time_start))

def extract_elements_from_list_by_index(targetList: list(), targetIndexList: List[int]):
    extractedList = []
    remainedList = []
    for i in range(len(targetList)):
        try:
            targetIndexList.index(i)
            extractedList.append(targetList[i])
        except ValueError:
            remainedList.append(targetList[i])
    return extractedList, remainedList

def edit_component_edges_by_germline(vertexList, edgeList, distThreshold, vGeneName, jGeneName):
    remainedVertexList = vertexList
    remainedEdgeList = edgeList
    resultEdgeList = []

    vjGeneName = vGeneName + '|' + jGeneName

    # distance from V genes
    distFromVgeneList = []
    for remainedVertex in remainedVertexList:
        distFromVgeneList.append(int(remainedVertex['distance_from_V_gene']))

    minDistFromVgene = min(distFromVgeneList)

    # change minDistFromVgene to str due to the limitation of read_file function
    tempToBeExtractedVertexIndexList = []
    for i, remainedVertex in enumerate(remainedVertexList):
        if remainedVertex['distance_from_V_gene'] == minDistFromVgene:
            tempToBeExtractedVertexIndexList.append(i)
            resultEdgeList.append({'from': vjGeneName, 'to':remainedVertex['name'], 'dist': minDistFromVgene+1})

    targetVertexList, remainedVertexList = extract_elements_from_list_by_index(remainedVertexList, tempToBeExtractedVertexIndexList)

    # Only two vertices and both shortest distance to the germline
    if len(targetVertexList) == 2 and len(remainedVertexList) == 0:
        resultEdgeList += remainedEdgeList
        return resultEdgeList

    while len(remainedVertexList) > 0:
        tempResultEdgeList = []

        targetVertexNameList = [targetVertex['name'] for targetVertex in targetVertexList]
        remainedVertexNameList = [remainedVertex['name'] for remainedVertex in remainedVertexList]

        tempToBeExtractedEdgeIndexList = []

        for i, edge in enumerate(remainedEdgeList):
            try:
                targetVertexNameList.index(edge['from'])
                tempToBeExtractedEdgeIndexList.append(i)
                continue
            except ValueError:
                pass

            try:
                targetVertexNameList.index(edge['to'])
                tempToBeExtractedEdgeIndexList.append(i)
                continue
            except ValueError:
                pass

        targetEdgeList, remainedEdgeList = extract_elements_from_list_by_index(remainedEdgeList, tempToBeExtractedEdgeIndexList)

        # eliminate closed loops
        tempToBeExtractedEdgeIndexList = []
        for i, targetEdge in enumerate(targetEdgeList):
            try:
                targetVertexNameList.index(targetEdge['from'])
                targetVertexNameList.index(targetEdge['to'])
                if targetEdge['dist'] <= distThreshold:
                    tempResultEdgeList.append(targetEdge)
            except ValueError:
                tempToBeExtractedEdgeIndexList.append(i)

        targetEdgeList, dummyEdgeList = extract_elements_from_list_by_index(targetList= targetEdgeList, targetIndexList= tempToBeExtractedEdgeIndexList)

        tempToBeExtractedVertexIndexList = []
        for i, targetEdge in enumerate(targetEdgeList):
            if targetEdge['dist'] <= distThreshold:
                tempResultEdgeList.append(targetEdge)
                try:
                    tempToBeExtractedVertexIndexList.append(remainedVertexNameList.index(targetEdge['from']))
                except ValueError:
                    tempToBeExtractedVertexIndexList.append(remainedVertexNameList.index(targetEdge['to']))

        tempToBeExtractedVertexIndexList = list(set(tempToBeExtractedVertexIndexList))

        if len(tempToBeExtractedVertexIndexList) > 0:
            targetVertexList, remainedVertexList = extract_elements_from_list_by_index(targetList= remainedVertexList, targetIndexList= tempToBeExtractedVertexIndexList)
            resultEdgeList += tempResultEdgeList
            continue

        else:
            targetEdgeDistList = [targetEdge['dist'] for targetEdge in targetEdgeList]
            tempToBeExtractedVertexIndexList = []
            minDist = min(targetEdgeDistList)
            for targetEdge in targetEdgeList:
                if targetEdge['dist'] == minDist:
                    tempResultEdgeList.append(targetEdge)
                    try:
                        tempToBeExtractedVertexIndexList.append(remainedVertexNameList.index(targetEdge['from']))
                    except ValueError:
                        tempToBeExtractedVertexIndexList.append(remainedVertexNameList.index(targetEdge['to']))

            tempToBeExtractedVertexIndexList = list(set(tempToBeExtractedVertexIndexList))

            targetVertexList, remainedVertexList = extract_elements_from_list_by_index(targetList= remainedVertexList, targetIndexList= tempToBeExtractedVertexIndexList)

            resultEdgeList += tempResultEdgeList

    for remainedEdge in remainedEdgeList:
        if remainedEdge['dist'] <= distThreshold:
            resultEdgeList.append(remainedEdge)

    return(resultEdgeList)

def edit_component_vertices_and_edges_by_germline(vertexHeader, vertexList, vjGeneSize, distThreshold):
    compVertexSet = {}
    compEdgeSet = {}

    vjGeneSet = {}

    resultVertexList = [x for x in vertexList]
    resultEdgeList = []

    for vertex in vertexList:
        compIndex = vertex['component_ID']
        vGene = vertex['V_gene'].split('*')[0]
        jGene = vertex['J_gene'].split('*')[0]
        vjGeneName = vGene + '|' + jGene

        try:
            prevVertexList = [x for x in compVertexSet[compIndex]]
            for prevVertex in prevVertexList:
                tempEdge = {'from': prevVertex['name'], 'to': vertex['name'],
                            'dist': Helper.levenshtein_distance(prevVertex['full_NT'], vertex['full_NT'])}
                compEdgeSet[compIndex].append(tempEdge)
            compVertexSet[compIndex].append(vertex)
        except KeyError:
            compVertexSet[compIndex] = [vertex]
            compEdgeSet[compIndex] = []

        try:
            vjGeneSet[vjGeneName] += 1
        except KeyError:
            vjGeneSet[vjGeneName] = 1

    for vjGene in vjGeneSet:
        tempVJGeneVertex = {}
        for colName in vertexHeader:
            tempVJGeneVertex[colName] = 'N/A'

        tempVJGeneVertex['name'] = vjGene
        tempVJGeneVertex['frequency'] = vjGeneSize

        resultVertexList.append(tempVJGeneVertex)

    for compIndex in compVertexSet:
        compVGeneName = compVertexSet[compIndex][0]['V_gene'].split('*')[0]
        compJGeneName = compVertexSet[compIndex][0]['J_gene'].split('*')[0]
        resultEdgeList += edit_component_edges_by_germline(vertexList= compVertexSet[compIndex],
                                                           edgeList= compEdgeSet[compIndex],
                                                           distThreshold= distThreshold,
                                                           vGeneName= compVGeneName,
                                                           jGeneName= compJGeneName)
    return resultVertexList, resultEdgeList


SAMPLE_NAME = 'SU0184'
print(SAMPLE_NAME)

DIR_CORRECTION = '/home/Human/data/AD/' + SAMPLE_NAME + '/2_error_correction'
DIR_ANNOTATION = '/home/Human/data/AD/' + SAMPLE_NAME + '/3_annotation'
DIR_FUNCTIONALITY = '/home/Human/data/AD/' + SAMPLE_NAME + '/4_functionality'
DIR_SAMPLING = '/home/Human/data/AD/' + SAMPLE_NAME + '/4-1_sampling'
DIR_NETWORK = '/home/Human/data/AD/' + SAMPLE_NAME + '/5_network'
IGHC_REFERENCE_FILE = '/home/Human/reference/IGHC_CH1.csv'

# Helper from the SeqUtil module
seqHelper = Helper()

if not os.path.exists(DIR_CORRECTION):
    os.makedirs(DIR_CORRECTION)

if not os.path.exists(DIR_ANNOTATION):
    os.makedirs(DIR_ANNOTATION)

if not os.path.exists(DIR_FUNCTIONALITY):
    os.makedirs(DIR_FUNCTIONALITY)

if not os.path.exists(DIR_NETWORK):
    os.makedirs(DIR_NETWORK)

if not os.path.exists(DIR_SAMPLING):
    os.makedirs(DIR_SAMPLING)

# Error correction
if False:
    # chain = ['VH']
    chainList = ['IPM_TH1240000']
    for chain in chainList:
        in_file = os.path.join(DIR_EXTRACTION, chain+'_p1_d1.csv')
        out_file = os.path.join(DIR_CORRECTION, chain+'_p1_d1_e1.csv')
        correct_error(file_in=in_file, file_out=out_file, evalue=0.001, col_sample=3, col_num=1, col_seq='full_NT')

# C gene annotation
if False:
    chRefHeder, chRefList = read_file(in_file=IGHC_REFERENCE_FILE)
    sampleList = [SAMPLE_NAME + '_d1_e1']

    isotypeList = ['M', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2']

    for sample in sampleList:
        targetSeqList = []
        readFileName = sample + '.csv'
        targetSeqHeader, targetSeqList = read_file(os.path.join(DIR_CORRECTION, readFileName))

        blastResultList = seqHelper.blast_run(query_list=[x['full_NT'] for x in targetSeqList], db_list=[x['full_NT'] for x in chRefList], m_type='NT', num_alignments=1, num_threads=3, blast_dir=DIR_ANNOTATION)

        annoWhole = 0
        annoSuccess = 0

        newSeqSet = {}

        for i, blastResult in enumerate(blastResultList):
            annoWhole += 1
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

                sbjctSeq = chRefList[int(alignment['sbjct_index'])]['full_NT']
                sbjctName = chRefList[int(alignment['sbjct_index'])]['name']
                querySeq = targetSeqList[i]['full_NT']
                queryRC = targetSeqList[i]['readcount']

                newSeq = querySeq[:alignment['query_start'] - alignment['sbjct_start']]
                newSeqIsotype = sbjctName.split('*')[0][3:]
                newSeqRC = queryRC

                try:
                    isotypeIndex = isotypeList.index(newSeqIsotype)
                    try:
                        newSeqSet[newSeq][isotypeIndex] += newSeqRC
                    except KeyError:
                        newSeqSet[newSeq] = [0 for i in range(len(isotypeList))]
                        newSeqSet[newSeq][isotypeIndex] += newSeqRC
                except ValueError:
                    pass
        with open(os.path.join(DIR_ANNOTATION, sample + '_a1.csv'), 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['full_NT', 'M', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2'])
            for newSeq in newSeqSet:
                writer.writerow([newSeq] + newSeqSet[newSeq])


# Annotation
if False:
    sampleList = [SAMPLE_NAME]

    annotationList = ['CDR1_NT', 'CDR2_NT', 'CDR3_NT', 'CDR1_AA', 'CDR2_AA', 'CDR3_AA',
                      'hit', 'V_gene', 'D_gene', 'J_gene', 'FR1_from', 'FR1_to']

    for sample in sampleList:
        targetSeqList = []
        rowList = []

        readFileName = sample + '_d1_e1_a1.csv'
        writeFileName = sample + 'd1_e1_a2.csv'

        seqHeader, seqList = read_file(os.path.join(DIR_ANNOTATION, readFileName))

        fullNTIndex = seqHeader.index('full_NT')

        for seq in seqList:
            targetSeqList.append(seq['full_NT'])

        resultList = seqHelper.annotate_regions(sequence_list=targetSeqList,
                                                chain_type=ChainType.HUMAN_HEAVY,
                                                file_blast=os.path.join(DIR_ANNOTATION, 'temp_igblast', sample + '_igblasted.txt'),
                                                annotation_list=annotationList)

        for seq, result in zip(seqList, resultList):
            result.update(seq)

        resultHeader = seqHeader[:fullNTIndex+1] + annotationList + seqHeader[fullNTIndex+1:]

        write_file_uncond(out_file=os.path.join(DIR_ANNOTATION, writeFileName), header=resultHeader, data=resultList)

# sequence trimming
if False:
    seqHeader, seqList = read_file(in_file=os.path.join(DIR_ANNOTATION, SAMPLE_NAME + '_d1_e1_a2.csv'))
    newSeqSet = {}
    isotypeList = ['M','G1','G2','G3','G4','A1','A2']
    for seq in seqList:
        if seq['FR1_from'] != 'N/A':
            newSeq = seq['full_NT'][seq['FR1_to']-int(seq['FR1_to']/3)*3:]
            try:
                for isotype in isotypeList:
                    newSeqSet[newSeq][isotype] += seq[isotype]
            except KeyError:
                seq['full_NT'] = newSeq
                seq['full_AA'] = seqHelper.translate(newSeq)
                newSeqSet[newSeq] = seq
    newSeqList = [newSeqSet[newSeq] for newSeq in newSeqSet]
    write_file_uncond(out_file = os.path.join(DIR_ANNOTATION, SAMPLE_NAME + '_d1_e1_a3.csv'), data = newSeqList, header = ['full_NT', 'full_AA', 'CDR1_NT', 'CDR2_NT', 'CDR3_NT', 'CDR1_AA', 'CDR2_AA', 'CDR3_AA', 'hit', 'V_gene', 'D_gene', 'J_gene', 'FR1_from', 'FR1_to', 'M', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2'])


# Rule out non-functional rows
if False:
    seqHelper = Helper()
    allIsotypeList = ['M', 'D', 'E', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2']
    for i in range(8):
        sampleIndex = i + 1
        readFileName = 'donor_' + str(sampleIndex) + '_prev_annotated.csv'
        writeFileName = 'donor_' + str(sampleIndex) + '_a1_f2.csv'

        seqHeader, seqList = read_file(in_file=os.path.join(DIR_ANNOTATION,readFileName))
        newSeqSet = {}
        for seq in seqList:
            fullNT = seq['full_NT']
            hit = seq['hit']
            vjFrame = seq['VJ_frame']
            vGene = seq['V_gene']
            dGene = seq['D_gene']
            jGene = seq['J_gene']
            fr1To = seq['FR1_to']
            isotype = seq['isotype']
            cdr1NT = seq['CDR1_NT']
            cdr1AA = seq['CDR1_AA']
            cdr2NT = seq['CDR2_NT']
            cdr2AA = seq['CDR2_AA']
            cdr3NT = seq['CDR3_NT']
            cdr3AA = seq['CDR3_AA']
            distFromVGene = seq['distance_from_V_gene']
            tpList = [seq['TP1'], seq['TP2'], seq['TP3'], seq['TP4'], seq['TP5']]

            if vGene == 'N/A' or jGene == 'N/A':
                continue

            if hit == False:
                continue

            if vjFrame == 'Out':
                continue

            if len(fullNT) < 250:
                continue

            if fr1To == 'N/A':
                continue

            if fr1To < 80:
                continue

            newFullNT = fullNT[fr1To-60:]
            newFullAA = seqHelper.translate(newFullNT)

            isotypeIndex = isotype[3:]

            try:
                for i, tp in enumerate(tpList):
                    newSeqSet[newFullNT][isotypeIndex][i] += tp
            except KeyError:
                newSeqSet[newFullNT] = {'full_NT':newFullNT, 'full_AA': seqHelper.translate(newFullNT),
                                        'CDR1_NT': cdr1NT, 'CDR1_AA': cdr1AA,
                                        'CDR2_NT': cdr2NT, 'CDR2_AA': cdr2AA,
                                        'CDR3_NT': cdr3NT, 'CDR3_AA': cdr3AA,
                                        'V_gene': vGene, 'D_gene': dGene, 'J_gene':jGene,
                                        'distance_from_V_gene': distFromVGene,
                                        'M': [0,0,0,0,0],
                                        'G1': [0,0,0,0,0],
                                        'G2': [0,0,0,0,0],
                                        'G3': [0,0,0,0,0],
                                        'G4': [0,0,0,0,0],
                                        'A1': [0,0,0,0,0],
                                        'A2': [0,0,0,0,0],
                                        'D': [0,0,0,0,0],
                                        'E': [0,0,0,0,0]}
                for i,tp in enumerate(tpList):
                    newSeqSet[newFullNT][isotypeIndex][i] += tp

        newSeqList = []

        for newSeq in newSeqSet:
            for allIsotype in allIsotypeList:
                tempSeq = newSeqSet[newSeq]
                tempFreq = [str(x) for x in tempSeq[allIsotype]]
                tempSeq[allIsotype] = '.'.join(tempFreq)

            newSeqList.append(tempSeq)

        newHeader = ['full_NT', 'full_AA', 'CDR1_NT', 'CDR1_AA', 'CDR2_NT', 'CDR2_AA', 'CDR3_NT', 'CDR3_AA',
                     'V_gene', 'D_gene', 'J_gene',
                     'distance_from_V_gene',
                     'M', 'D', 'E', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2']

        write_file_uncond(out_file = os.path.join(DIR_FUNCTIONALITY, writeFileName), data = newSeqList, header = newHeader)

# annotate additional info.
if False:
    preProcessingType = 'd1_e1_a1_f2'
    afterProcessingType = 'd1_e1_a2_f2'

    targetSeqList = []
    annotationList = ['distance_from_V_gene', 'aligned_length_to_V_gene', 'FR3_to']

    seqHeader, seqList = read_file(in_file = os.path.join(DIR_FUNCTIONALITY, SAMPLE_NAME + '_' + preProcessingType + '.csv'))

    targetSeqList = [seq['full_NT'] for seq in seqList]

    annotatedList = seqHelper.annotate_regions(chain_type=ChainType.HUMAN_HEAVY, sequence_list=targetSeqList, annotation_list=annotationList)

    for seq, annotated in zip(seqList, annotatedList):
        seq['FR3_to'] = annotated['FR3_to']
        seq['distance_from_V_gene'] = annotated['distance_from_V_gene']
        seq['aligned_length_to_V_gene'] = annotated['aligned_length_to_V_gene']

    seqHeader.insert(seqHeader.index('FR1_to'), 'distance_from_V_gene')
    seqHeader.insert(seqHeader.index('FR1_to'), 'aligned_length_to_V_gene')
    seqHeader.insert(seqHeader.index('FR1_to'), 'FR3_to')

    seqHeader.remove('FR1_to')

    write_file_uncond(out_file=os.path.join(DIR_FUNCTIONALITY, SAMPLE_NAME + '_' + afterProcessingType + '.csv'),
                      header=seqHeader, data=seqList)

# Sampling for stat.
if False:
    preProcessingType = 'd1_e1_a2_f2'
    targetFile = os.path.join(DIR_FUNCTIONALITY, '_'.join([SAMPLE_NAME, preProcessingType + '.csv']))

    samplingNum = 160743

    targetHeader, targetList = read_file(in_file=targetFile)

    writeFileName = os.path.basename(targetFile).split('.')[0] + '_s1.csv'

    tempSampleList = []

    for target in targetList:
        rc = int(target['readcount'])
        for i in range(rc):
            tempSampleList.append(target)
    
    samplingResultList = random.sample(tempSampleList, samplingNum)

    samplingResultSet = {}

    for samplingResult in samplingResultList:
        try:
            samplingResultSet[str(samplingResult)].append(samplingResult)
        except KeyError:
            samplingResultSet[str(samplingResult)] = [samplingResult]

    samplingResultList = []
    for samplingResult in samplingResultSet:
        tempResult = samplingResultSet[samplingResult][0]
        tempResult['readcount'] = len(samplingResultSet[samplingResult])
        samplingResultList.append(tempResult)

    write_file_uncond(out_file=os.path.join(DIR_SAMPLING, writeFileName), header=targetHeader, data=samplingResultList)

# Make vertex files
if False:
    preProcessingType = 'd1_e1_a2_f2_s1'
    vertexType = 'nv1'

    file_target = os.path.join(DIR_SAMPLING, '_'.join([SAMPLE_NAME, preProcessingType + '.csv']))
    file_output = os.path.join(DIR_NETWORK, '_'.join([SAMPLE_NAME, preProcessingType, vertexType + '.csv']))
    set_id(file_target, file_output, col_id='seq_id', col_seq='full_NT', id_prefix=SAMPLE_NAME + '-')

# Extract edges from vertex files
if False:
    sample = SAMPLE_NAME
    preProcessingType = 'd1_e1_a2_f2_s1'
    threshold = 8
    vertexType = 'nv1'
    edgeType = 'ne8cdr'
    regionType = 'CDR3_NT'

    file_target = os.path.join(DIR_NETWORK, '_'.join([sample, preProcessingType, vertexType])+'.csv')
    file_output = os.path.join(DIR_NETWORK, '_'.join([sample, preProcessingType, vertexType, edgeType])+'.csv')
    extract_edge_all_levenstein_under_thres(file_in_nodes=file_target,
                                            file_out_edges=file_output,
                                            col_id='seq_id',
                                            col_seq=regionType,
                                            threshold=threshold,
                                            num_process=8)

# calculate minimum distance of the vertices for NN plot
if False:
    sample = SAMPLE_NAME
    preProcessingType = 'd1_e1_a2_f2_s1'
    vertexType = 'nv1'
    edgeType = 'ne8cdr'

    vertexColNameSeqID = 'seq_id'
    vertexColNameSeq = 'CDR3_NT'

    edgeColNameFrom = 'from'
    edgeColNameTo = 'to'
    edgeColNameDist = 'dist'
    edgeColNameSeqLength = 'sequence_length'

    minDistSet = {}
    seqLengthSet = {}

    targetVertexFileName = '_'.join([sample, preProcessingType, vertexType]) + '.csv'
    targetEdgeFileName = '_'.join([sample, preProcessingType, vertexType, edgeType]) + '.csv'

    writeMinDistFileName = '_'.join([sample, preProcessingType, vertexType, edgeType, 'NN_plot']) + '.csv'

    targetVertexFile = os.path.join(DIR_NETWORK, targetVertexFileName)
    targetEdgeFile = os.path.join(DIR_NETWORK, targetEdgeFileName)

    writeMinDistFile = os.path.join(DIR_NETWORK, writeMinDistFileName)

    vertexHeader, vertexList = read_file(in_file=targetVertexFile)
    print('vertex file reading done!')

    # edgeHeader, edgeList = read_file(in_file=targetEdgeFile)
    edgeList = []
    with open(targetEdgeFile, 'r') as f:
        reader = csv.reader(f)
        edgeHeader = next(reader)
        for row in reader:
            tempSet = {}
            tempSet['from'] = row[0]
            tempSet['to'] = row[1]
            tempSet['dist'] = int(row[2])
            if tempSet['dist'] > 0:
                edgeList.append(tempSet)
    print('edge file reading done!')

    print(sample + ' min dist extraction started!')

    for i, vertex in enumerate(vertexList):
        seqID = vertex[vertexColNameSeqID]
        minDistSet[seqID] = []
        seqLengthSet[seqID] = len(vertex[vertexColNameSeq])

    for edge in edgeList:
        edgeFrom = edge[edgeColNameFrom]
        edgeTo = edge[edgeColNameTo]
        edgeDist = int(edge[edgeColNameDist])

        minDistSet[edgeFrom].append(edgeDist)
        minDistSet[edgeTo].append(edgeDist)

    for i, seqID in enumerate(minDistSet):
        if len(minDistSet[seqID]) == 0:
            minDistSet[seqID] = ['N/A', seqLengthSet[seqID]]
        else:
            # minDistSet[seqID] = min(minDistSet[seqID]) / seqLengthSet[seqID]
            minDistSet[seqID] = [min(minDistSet[seqID]), seqLengthSet[seqID]]

    with open(writeMinDistFile, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['seq_id', 'min_dist', 'seq_len'])
        for seqID in minDistSet:
            writer.writerow([seqID] + minDistSet[seqID])

# edit edge file to utilize NN plot results & gene usage information
if False:
    sample = SAMPLE_NAME
    preProcessingType = 'd1_e1_a2_f2_s1'
    vertexType = 'nv1'
    edgeType = 'ne8cdr'

    sampleList = ['SU0040', 'SU0101', 'SU0104', 'SU0185', 'SU0226', 'SU0250', 'SU0304', 'SU0428', 'SU0460']
    CDR3LDThresholdList = [0.09205, 0.08395, 0.06146, 0.06958, 0.06593, 0.06797, 0.07850, 0.06687, 0.04997]

    sampleList = ['SU0313', 'SU0292', 'SU0202', 'SU0184']
    CDR3LDThresholdList = [0.050906, 0.06184, 0.091778, 0.079944]

    CDR3LDThreshold = CDR3LDThresholdList[sampleList.index(SAMPLE_NAME)]

    lengthDiffThreshold = 3
    writeEdgeType = 'ne8a'

    print(sample + ' edge edit was started!')
    vertexFileName = '_'.join([sample, preProcessingType, vertexType]) + '.csv'
    edgeFileName = '_'.join([sample, preProcessingType, vertexType, edgeType]) + '.csv'
    vertexFile = os.path.join(DIR_NETWORK, vertexFileName)
    edgeFile = os.path.join(DIR_NETWORK, edgeFileName)

    vertexHeader, vertexList = read_file(vertexFile)

    edgeList = []
    with open(edgeFile, 'r') as f:
        reader =csv.reader(f)
        edgeHeader = next(reader)
        for row in reader:
            edgeList.append({'from':row[0], 'to':row[1], 'dist':int(row[2])})

    print('Vertex & Edge file reading was done!')

    newEdgeList = []

    vertexReferenceSet = {}

    for vertex in vertexList:
        seqName = vertex['seq_id']
        vertexReferenceSet[seqName] = vertex

    for i, edge in enumerate(edgeList):
        fromVertexID = edge['from']
        toVertexID = edge['to']
        distance = edge['dist']

        fromVertexSeqLength = len(vertexReferenceSet[fromVertexID]['CDR3_NT'])
        toVertexSeqLength = len(vertexReferenceSet[toVertexID]['CDR3_NT'])
        averageSeqLength = np.mean([fromVertexSeqLength, toVertexSeqLength])

        if float(distance) / averageSeqLength > CDR3LDThreshold:
            continue

        fromVertex = vertexReferenceSet[fromVertexID]
        toVertex = vertexReferenceSet[toVertexID]

        fromVertexVGene = fromVertex['V_gene']
        fromVertexJGene = fromVertex['J_gene']

        toVertexVGene = toVertex['V_gene']
        toVertexJGene = toVertex['J_gene']

        if fromVertexVGene != toVertexVGene:
            continue

        if fromVertexJGene != toVertexJGene:
            continue

        fromVertexFullNT = fromVertex['full_NT']
        toVertexFullNT = toVertex['full_NT']

        if abs(len(fromVertexFullNT) - len(toVertexFullNT)) > lengthDiffThreshold:
            continue

        newEdgeList.append(edge)

    writeEdgeFileName = '_'.join([sample, preProcessingType, vertexType, writeEdgeType]) + '.csv'
    writeEdgeFile = os.path.join(DIR_NETWORK, writeEdgeFileName)

    write_file_uncond(out_file=writeEdgeFile, header=edgeHeader, data=newEdgeList)
