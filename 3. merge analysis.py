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
import numpy as np

from Common.FileUtil import read_file, write_file, write_file_uncond

from Common.SeqUtil_test import Helper
from typing import List, Dict, Tuple
import csv

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

    write_file_uncond(file_out, header, data)

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


workDir = '/home/Human/data/AD/files_for_stat/merge_analysis/pairwise_cdr3aa_merge_analysis/pair_data'

pairList = ['11&12','11&13','11&21','11&22','11&23','11&31','11&32','11&33',
            '12&13','12&21','12&22','12&23','12&31','12&32','12&33',
            '13&21','13&22','13&23','13&31','13&32','13&33',
            '21&22','21&23','21&31','21&32','21&33',
            '22&23','22&31','22&32','22&33',
            '23&31','23&32','23&33',
            '31&32','31&33',
            '32&33']

isotypeList = ['M', 'D', 'E', 'G1', 'G2', 'G3', 'G4', 'A1', 'A2']

# Merge data
if False:
    sampleSet = {'11':'SU0304', '12':'SU0101', '13':'SU0250',
                 '21':'SU0226', '22':'SU0460', '23':'SU0428',
                 '31':'SU0040', '32':'SU0104', '33':'SU0185',
                 '41':'SU0184', '42':'SU0202', '43':'SU0292',
                 '52':'SU0313'}

    sampleList = list(sampleSet.keys())

    addList = ['41', '42', '43', '52']

    workDir = '/home/Human/data/AD'

    for i in range(len(sampleList)):
        refSample = sampleList[i]
        querySampleList = sampleList[i+1:]
        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue

            print(refSample, querySample)

            targetDir1 = os.path.join(workDir, sampleSet[refSample], '5_network')
            targetDir2 = os.path.join(workDir, sampleSet[querySample], '5_network')
            if refSample == '41' or refSample == '42' or refSample == '43' or refSample =='52':
                targetFileName1 = sampleSet[refSample] + '_d1_e1_a2_f2_s1_nv1a.csv'
            else:
                targetFileName1 = sampleSet[refSample] + '_d1_e1_a2_f1_s1_nv1a.csv'

            if querySample == '41' or querySample == '42' or querySample == '43' or querySample =='52':
                targetFileName2 = sampleSet[querySample] + '_d1_e1_a2_f2_s1_nv1a.csv'
            else:
                targetFileName2 = sampleSet[querySample] + '_d1_e1_a2_f1_s1_nv1a.csv'

            targetHeader1, targetList1 = read_file(os.path.join(targetDir1, targetFileName1))
            targetHeader2, targetList2 = read_file(os.path.join(targetDir2, targetFileName2))

            writeList = []

            for target in targetList1:
                target['sample'] = refSample
                writeList.append(target)

            for target in targetList2:
                target['sample'] = querySample
                writeList.append(target)

            writeSample = refSample + '&' + querySample

            writeDir = os.path.join('/home/Human/data/AD/files_for_stat/merge_analysis/pairwise_cdr3aa_merge_analysis/pair_data', writeSample)

            if os.path.exists(writeDir) == False:
                os.makedirs(writeDir)

            writeFileName = writeSample + '.csv'

            write_file_uncond(out_file=os.path.join(writeDir, writeFileName),
                              header=targetHeader1 + ['sample'],
                              data=writeList)

# Make vertex files
if False:
    sampleList = ['11', '12', '13',
                  '21', '22', '23',
                  '31', '32', '33',
                  '41', '42', '43',
                  '52']

    addList = ['41', '42', '43', '52']

    for i in range(len(sampleList)):
        refSample = sampleList[i]

        querySampleList = sampleList[i+1:]

        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue
            pair = refSample + '&' + querySample
            print(pair)
            SAMPLE_NAME = pair
            sampleDir = os.path.join(workDir, pair)
            vertexType = 'nv2'

            file_target = os.path.join(sampleDir, '_'.join([SAMPLE_NAME + '.csv']))
            file_output = os.path.join(sampleDir, '_'.join([SAMPLE_NAME, vertexType + '.csv']))
            set_id(file_target, file_output, col_id='seq_id', col_seq='full_NT', id_prefix=SAMPLE_NAME + '-')

# Extract edges from vertex files
if False:
    sampleList = ['11', '12', '13',
                  '21', '22', '23',
                  '31', '32', '33',
                  '41', '42', '43',
                  '52']

    addList = ['41', '42', '43', '52']

    for i in range(len(sampleList)):
        refSample = sampleList[i]

        querySampleList = sampleList[i + 1:]

        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue

            pair = refSample + '&' + querySample

            print(pair)

            sampleDir = os.path.join(workDir, pair)

            sample = pair
            threshold = 6
            vertexType = 'nv2'
            edgeType = 'ne6cdr3aa'
            regionType = 'CDR3_AA'

            file_target = os.path.join(sampleDir, '_'.join([sample, vertexType])+'.csv')
            file_output = os.path.join(sampleDir, '_'.join([sample, vertexType, edgeType])+'.csv')
            extract_edge_all_levenstein_under_thres(file_in_nodes=file_target,
                                                    file_out_edges=file_output,
                                                    col_id='seq_id',
                                                    col_seq=regionType,
                                                    threshold=threshold,
                                                    num_process=30)

# calculate minimum distance of the vertices for NN plot
if False:
    sampleList = ['11', '12', '13',
                  '21', '22', '23',
                  '31', '32', '33',
                  '41', '42', '43',
                  '52']

    addList = ['41', '42', '43', '52']

    for i in range(len(sampleList)):
        refSample = sampleList[i]

        querySampleList = sampleList[i + 1:]

        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue

            pair = refSample + '&' + querySample

            print(pair)

            sampleDir = os.path.join(workDir, pair)
            sample = pair
            vertexType = 'nv2'
            edgeType = 'ne6cdr3aa'

            vertexColNameSeqID = 'seq_id'
            vertexColNameSeq = 'CDR3_AA'

            edgeColNameFrom = 'from'
            edgeColNameTo = 'to'
            edgeColNameDist = 'dist'
            edgeColNameSeqLength = 'sequence_length'

            minDistSet = {}
            seqLengthSet = {}

            targetVertexFileName = '_'.join([sample, vertexType]) + '.csv'
            targetEdgeFileName = '_'.join([sample, vertexType, edgeType]) + '.csv'

            writeMinDistFileName = '_'.join([sample, vertexType, edgeType, 'NN_plot']) + '.csv'

            targetVertexFile = os.path.join(sampleDir, targetVertexFileName)
            targetEdgeFile = os.path.join(sampleDir, targetEdgeFileName)

            writeMinDistFile = os.path.join(sampleDir, writeMinDistFileName)

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

# convert to net file
if False:
    sampleList = ['11', '12', '13',
                  '21', '22', '23',
                  '31', '32', '33',
                  '41', '42', '43',
                  '52']

    addList = ['41', '42', '43', '52']

    for i in range(len(sampleList)):
        refSample = sampleList[i]

        querySampleList = sampleList[i + 1:]

        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue

            pair = refSample + '&' + querySample

            print(pair)
            sampleDir = os.path.join(workDir, pair)
            sample = pair
            vertexType = 'nv1'
            writeVertexType = 'nv2'

            readFileName = '_'.join([sample, vertexType + '.csv'])
            writeFileName = '_'.join([sample, writeVertexType + '.csv'])

            targetHeader, targetList = read_file(in_file=os.path.join(sampleDir, readFileName))

            fullNTSet = {}
            for target in targetList:
                tempKey = str(target['sample']) + '|' + target['full_NT']
                try:
                    fullNTSet[tempKey][target['isotype']] += int(target['readcount'])
                except KeyError:
                    fullNTSet[tempKey] = target
                    for isotype in isotypeList:
                        fullNTSet[tempKey][isotype] = 0
                    fullNTSet[tempKey][target['isotype']] += int(target['readcount'])

            writeList = list(fullNTSet.values())

            targetHeader.remove('readcount')
            targetHeader.remove('isotype')
            targetHeader += isotypeList

            write_file_uncond(out_file=os.path.join(sampleDir, writeFileName),
                              header=targetHeader,
                              data=writeList)


def cluster_components_by_cdr3aa_length_and_vjGene(compSet):
    tempSet = {}
    for compID in compSet:
        tempVertex = compSet[compID][0]
        cdr3aaLen = len(tempVertex['CDR3_AA'])
        vGene = tempVertex['V_gene'].split('*')[0]
        jGene = tempVertex['J_gene'].split('*')[0]
        tempKey = '|'.join([str(cdr3aaLen), str(vGene), str(jGene)])
        try:
            tempSet[tempKey].append(compSet[compID])
        except KeyError:
            tempSet[tempKey] = [compSet[compID]]
    clusteredList= list(tempSet.values())
    return clusteredList

def sort_by_readcount_sum(compList, isotypeList):
    tempSet = {}
    for comp in compList:
        compRcSum = 0
        for vertex in comp:
            for isotype in isotypeList:
                compRcSum += int(vertex[isotype])
        try:
            tempSet[compRcSum].append(comp)
        except KeyError:
            tempSet[compRcSum] = [comp]

    rcSumKeyList = list(tempSet.keys())
    rcSumKeyList.sort(reverse=True)

    sortedCompList = []
    for rcSumKey in rcSumKeyList:
        sortedCompList += tempSet[rcSumKey]

    return sortedCompList

def get_averaged_cdr3aa_distance(comp1, comp2):
    tempHelper = Helper()
    distList = []
    for vertex1 in comp1:
        cdr3aa1 = vertex1['CDR3_AA']
        for vertex2 in comp2:
            cdr3aa2 = vertex2['CDR3_AA']
            cdr3aaDist = tempHelper.hamming_distance(cdr3aa1, cdr3aa2)
            distList.append(cdr3aaDist)
    averagedDist = np.mean(distList)
    return averagedDist

def get_min_cdr3aa_distance(comp1, comp2):
    tempHelper = Helper()
    distList = []
    for vertex1 in comp1:
        cdr3aa1 = vertex1['CDR3_AA']
        for vertex2 in comp2:
            cdr3aa2 = vertex2['CDR3_AA']
            cdr3aaDist = tempHelper.hamming_distance(cdr3aa1, cdr3aa2)
            distList.append(cdr3aaDist)
    minDist = np.min(distList)
    return minDist

def merge_components_by_averaged_cdr3aa_distance(compList, distThresh):
   if len(compList) == 0:
       print('Empty list for cdr3aa distance component merging.')
       raise RuntimeError

   elif len(compList) == 1:
       return compList

   else:
       resultCompList = []
       while len(compList) > 0:
           refComp = compList.pop()

           if len(compList) == 0:
               resultCompList.append(refComp)
               break

           else:
               removeIndexList = []
               for i, comp in enumerate(compList):
                   cdr3aaDist = get_averaged_cdr3aa_distance(refComp, comp)
                   if cdr3aaDist < distThresh:
                       refComp += comp
                       removeIndexList.insert(0, i)
               for removeIndex in removeIndexList:
                   del compList[removeIndex]

               resultCompList.append(refComp)

       return resultCompList

def merge_components_by_min_cdr3aa_distance(compList, distThresh):
   if len(compList) == 0:
       print('Empty list for cdr3aa distance component merging.')
       raise RuntimeError

   elif len(compList) == 1:
       return compList

   else:
       resultCompList = []
       while len(compList) > 0:
           refComp = compList.pop()

           if len(compList) == 0:
               resultCompList.append(refComp)
               break

           else:
               removeIndexList = []
               for i, comp in enumerate(compList):
                   cdr3aaDist = get_min_cdr3aa_distance(refComp, comp)
                   if cdr3aaDist < distThresh:
                       refComp += comp
                       removeIndexList.insert(0, i)
               for removeIndex in removeIndexList:
                   del compList[removeIndex]

               resultCompList.append(refComp)

       return resultCompList


def sort_components_by_its_size(compList):
    compSet = {}
    for comp in compList:
        compSize = len(comp)
        try:
            compSet[compSize].append(comp)
        except KeyError:
            compSet[compSize] = [comp]
    compSizeKeyList = list(compSet.keys())
    compSizeKeyList.sort(reverse=True)
    sortedCompList = []
    for compSizeKey in compSizeKeyList:
        sortedCompList += compSet[compSizeKey]
    return sortedCompList

# component merging approach
if False:
    distThreshList = [0.128281, 0.131505, 0.129808, 0.10875,
                      0.135156, 0.140741, 0.133333, 0.11111,
                      0.127734, 0.130862, 0.130862, 0.108621,

                      0.128143, 0.131375, 0.130571, 0.108714,
                      0.121797, 0.118889, 0.123704, 0.106852,
                      0.128125, 0.133448, 0.128621, 0.111724,

                      0.135156, 0.140741, 0.134712, 0.119904,
                      0.120312, 0.124643, 0.124643, 0.102411,
                      0.127734, 0.123448, 0.130862, 0.108621,

                      0.135391, 0.130547, 0.120859,
                      0.133333, 0.116481,
                      0.122596]

    distThreshList = [0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,

                      0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,

                      0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,
                      0.3, 0.3, 0.3, 0.3,

                      0.3, 0.3, 0.3,
                      0.3, 0.3,
                      0.3]

    pairList = ['11&41','11&42','11&43','11&52',
                '12&41','12&42','12&43','12&52',
                '13&41','13&42','13&43','13&52',

                '21&41','21&42','21&43','21&52',
                '22&41','22&42','22&43','22&52',
                '23&41','23&42','23&43','23&52',

                '31&41', '31&42', '31&43', '31&52',
                '32&41', '32&42', '32&43', '32&52',
                '33&41', '33&42', '33&43', '33&52',

                '41&42','41&43','41&52',
                '42&43','42&52',
                '43&52']

    for pair, distThresh in zip(pairList, distThreshList):
        print(pair)
        sampleDir = os.path.join(workDir, pair)
        sample = pair
        vertexType = 'nv2'
        edgeType = 'ne6cdr3aa'

        writeVertexType = 'nv2b'

        sampleFileName = '_'.join([sample, vertexType]) + '.csv'
        writeFileName = '_'.join([sample, writeVertexType]) + '.csv'

        mergeHeader, mergeList = read_file(in_file=os.path.join(sampleDir, sampleFileName))

        compSet = {}

        for merge in mergeList:
            compID = str(merge['sample']) + '|' + str(merge['component_ID'])
            try:
                compSet[compID].append(merge)
            except KeyError:
                compSet[compID] = [merge]

        clusteredCompList = cluster_components_by_cdr3aa_length_and_vjGene(compSet)

        mergedCompList = []
        for clusteredComps in clusteredCompList:
            sortedCompList = sort_by_readcount_sum(compList=clusteredComps, isotypeList=isotypeList)

            # tempMergedCompList = merge_components_by_averaged_cdr3aa_distance(sortedCompList, distThresh)
            tempMergedCompList = merge_components_by_min_cdr3aa_distance(sortedCompList, distThresh)
            mergedCompList += tempMergedCompList

        sortedMergedCompList = sort_components_by_its_size(mergedCompList)

        for i, sortedMergedComp in enumerate(sortedMergedCompList):
            for vertex in sortedMergedComp:
                vertex['component_ID_merged'] = i + 1

        resultVertexList = []
        for sortedMergedComp in sortedMergedCompList:
            resultVertexList += sortedMergedComp

        write_file_uncond(out_file=os.path.join(sampleDir, writeFileName),
                          header=mergeHeader + ['component_ID_merged'],
                          data=resultVertexList)


# Extract shared components
if False:
    targetSampleList = ['11', '12', '13',
                  '21', '22', '23',
                  '31', '32', '33',
                  '41', '42', '43',
                  '52']

    addList = ['41', '42', '43', '52']

    for i in range(len(targetSampleList)):
        refSample = targetSampleList[i]

        querySampleList = targetSampleList[i + 1:]

        for querySample in querySampleList:
            try:
                addList.index(refSample)
            except ValueError:
                try:
                    addList.index(querySample)
                except ValueError:
                    continue

            pair = refSample + '&' + querySample

            print(pair)
            sample = pair
            sampleDir = os.path.join(workDir, pair)
            writeDir = '/home/Human/data/AD/files_for_stat/merge_analysis/pairwise_cdr3aa_merge_analysis/shared_component'
            readFileName = sample + '_nv2b.csv'
            writeFileName = sample + '.csv'

            targetHeader, targetList = read_file(in_file=os.path.join(sampleDir, readFileName))

            compSet = {}
            for target in targetList:
                try:
                    compSet[target['component_ID_merged']].append(target)
                except KeyError:
                    compSet[target['component_ID_merged']] = [target]

            resultList = []
            for comp in compSet:
                sampleList = []
                vertexList = compSet[comp]
                for vertex in vertexList:
                    sampleList.append(vertex['sample'])
                sampleList = list(set(sampleList))
                if len(sampleList) > 1:
                    resultList += vertexList

            write_file_uncond(out_file=os.path.join(writeDir, writeFileName),
                              header=targetHeader,
                              data=resultList)

def get_min_dist_between_comps(comp1, comp2):
    minDist = float('inf')
    tempHelper = Helper()
    for vertex1 in comp1:
        for vertex2 in comp2:
            dist = tempHelper.levenshtein_distance(vertex1['full_AA'], vertex2['full_AA'])
            if dist < minDist:
                minDist = dist
    return minDist

def get_average_dist_between_comps(comp1, comp2):
    tempHelper = Helper()
    distList = []
    for vertex1 in comp1:
        for vertex2 in comp2:
            dist = tempHelper.levenshtein_distance(vertex1['full_AA'], vertex2['full_AA'])
            distList.append(dist)
    averageDist = np.mean(distList)
    return averageDist



