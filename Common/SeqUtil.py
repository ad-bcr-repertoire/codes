import copy
import logging
import os
import re
import tempfile
from multiprocessing import Process
from pathlib import Path

import editdistance
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACProtein, ExtendedIUPACProtein
from Bio.Data import CodonTable
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import List, Dict, Tuple

from Common.Enum import ChainType
from Common.FileUtil import detach_format, TableData, TableHeader
from Common.ProcessUtil import run_igblast

# New import for blast_run
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIStandalone

_protein_alphabet = IUPACProtein()
_protein_letters = str(_protein_alphabet.letters)


def is_nucleotide_or_protein(sequence: str):
    """
    return TRUE if the input sequence is a dna sequence,
    FALSE if it is a protein sequence,
    None otherwise.
    :param sequence: 
    :return: 
    """
    if all(ch.upper() in 'AGTC-. ' for ch in sequence):
        return True
    elif all(ch.upper() in _protein_letters + '-. X' for ch in sequence):
        return False
    else:
        return None


class PTMFinder(object):
    """Post Translation Modification sites finder.
    """

    def __init__(self):
        """"""
        self.re_deamidation = re.compile('(NS|NG|NH)')
        self.re_aspartate_isomerization = re.compile('(DG|DS|DT|DD|DA)')
        self.re_N_glycosylation = re.compile('(N[^P][ST])')
        self.re_cleavage = re.compile('(DP|DQ|NS)')
        self.re_oxidation = re.compile('([MW])')
        self.re_free_cysteine = re.compile('(C)')
        self.re_sulfation = re.compile('(Y)')
        self.re_methylation = re.compile('([KR])')

    def find_diamidation(self, sequence) -> List:
        """
        Find deamidation sites from the given sequence.
        
        * Deamidation Sites    
        Asparagine residues as potential deamidation sites were analysed in the context of
        both the amino and carboxy adjacent amino acids using a method based on
        Robinson & Robinson, PNAS (2001) 98,3, p944-949.
        No potential deamidation sites with T1/2 <25 days were identified within the VL.
        
        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_deamidation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_deamidation.finditer(str(sequence))]
        else:
            return []

    def find_aspartate_isomerization(self, sequence) -> List:
        """
        Find aspartate isomerization sites from the given sequence.

        * Aspartate Isomerization Sites    
        Aspartate isomerization sites were predicted by analysing the sequence
        for known isomerization motifs: DG, DS, DT, DD or DA.
        Two potential isomerization sites were identified in the VL sequence.
        Aspartate 91 and Aspartate 95A, both of which are located in CDR3.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_aspartate_isomerization.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_aspartate_isomerization.finditer(str(sequence))]
        else:
            return []

    def find_n_glycosylation(self, sequence) -> List:
        """
        Find N glycosylation sites from the given sequence.

        * N Glycosylation Sites    
        Sequences were analysed based on the consensus N-linked glycosylation motif:
        - N-X-S/T - (where X can be any amino acid except Proline)
        No potential N-linked glycosylation sites were identified within the VL domain.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_N_glycosylation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_N_glycosylation.finditer(str(sequence))]
        else:
            return []

    def find_cleavage(self, sequence) -> List:
        """
        Find cleavage sites from the given sequence.
        
        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_cleavage.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_cleavage.finditer(str(sequence))]
        else:
            return []

    def find_oxidation(self, sequence) -> List:
        """
        Find oxidation sites from the given sequence.

        * Oxidation Sites    
        Structural models were prepared and methionine and tryptophan residues were identified and assessed
        to determine whether they are likely to be surface exposed (and therefore candidates for oxidation).

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_oxidation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_oxidation.finditer(str(sequence))]
        else:
            return []

    def find_free_cysteine(self, sequence) -> List:
        """
        Find free cysteine sites from the given sequence.

        * Free Cysteine Sites Sequences were analysed to identify unpaired cysteine residues. VL domains normally 
        contain two cysteines which form a disulphide bond in the folded molecule. Additional cysteines would be 
        expected to be detrimental to folding and potentially cause issues such as aggregation. No unpaired cysteines 
        were identified within the VL domain. 

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_free_cysteine.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_free_cysteine.finditer(str(sequence))]
        else:
            return []

    def find_sulfation(self, sequence) -> List:
        """
        Find sulfation sites from the given sequence.

        * Sulfation Sites Tyrosine sulfation sites were predicted using the method of Monigatti F. et al., 
        Bioinformatics 18:769-770 (2002). No sulfated tyrosine sites were predicted 

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_sulfation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_sulfation.finditer(str(sequence))]
        else:
            return []

    def find_methylation(self, sequence) -> List:
        """
        Find methylation sites from the given sequence.

        * Methylation Sites    
        Prediction of methylation of lysine and arginine residues was performed using the method of
        Chen et al., Nucleic Acids Research  34 (Issue suppl 2): 249-253 (2006)
        One methylated lysine site was predicted – Lysine 66 (located in Framework 3)
        PSRFSGS - K - SGSTHTL
        No methylated arginine sites were predicted.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in self.re_methylation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in self.re_methylation.finditer(str(sequence))]
        else:
            return []

    @staticmethod
    def find_custom(sequence, re_compiled) -> List:
        """
        Find custom subsequence  sites from the given sequence.
        
        :param sequence: Target amino acid sequence to find.
        :param re_compiled: 
        :return: List of match objects that refer the sites found. 
        """

        if type(sequence) is str:
            return [_ for _ in re_compiled.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (ExtendedIUPACProtein, IUPACProtein):
            return [_ for _ in re_compiled.finditer(str(sequence))]
        else:
            return []


class Helper(object):
    def __init__(self):
        """"""
        self.reverse_sequence = {
            'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G',
            'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
            'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
            'N': 'N', ' ': ' ', '.': '.'
        }
        self.standard_table = CodonTable.unambiguous_dna_by_id[1]
        self.standard_codon_map = copy.deepcopy(self.standard_table.forward_table)
        for triplet in self.standard_table.stop_codons:
            self.standard_codon_map[triplet] = '*'
        self.standard_codon_map['   '] = ' '
        self.standard_codon_map['...'] = '.'
        self.dna_symbol_map = {
            'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U',
            'W': '[AT]', 'S': '[CG]', 'M': '[AC]', 'K': '[GT]', 'R': '[AG]', 'Y': '[CT]',
            'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
            'N': '[ACGT]'
        }

        # amber codon TAG->Q
        self.amber_codon_map = copy.deepcopy(self.standard_table.forward_table)
        for triplet in self.standard_table.stop_codons:
            self.amber_codon_map[triplet] = '*'
        self.amber_codon_map['TAG'] = 'Q'

    def to_reverse_complement(self, sequence: str) -> str:
        """"""
        return ''.join(self.reverse_sequence[s] for s in sequence)[::-1]

    def to_regex_string_from_dna(self, sequence: str, is_forward: bool = True) -> str:
        """
        Symbol Table
        A	Adenine	A
        C	Cytosine		C
        G	Guanine			G
        T	Thymine				T
        U	Uracil				U
        W	Weak	A			T
        S	Strong		C	G
        M	aMino	A	C
        K	Keto			G	T
        R	puRine	A		G	
        Y	pYrimidine		C		T
        B	not A (B comes after A)		C	G	T	
        D	not C (D comes after C)	A		G	T
        H	not G (H comes after G)	A	C		T
        V	not T (V comes after T and U)	A	C	G	
        N	any Nucleotide (not a gap)	A	C	G	T	
        :param sequence: 
        :return: 
        """

        if is_forward:
            return ''.join([self.dna_symbol_map[c] for c in sequence.upper()])
        else:
            return ''.join([self.dna_symbol_map[c] for c in self.to_reverse_complement(sequence.upper())])

    @staticmethod
    def hamming_distance(s1: str, s2: str) -> int:
        """Calculate the Hamming distance between two strings"""
        if len(s1) != len(s2):
            return max(len(s1), len(s2))

        return sum((c1 != c2 and c1 != '.' and c2 != '.') for c1, c2 in zip(s1, s2))

    @staticmethod
    def get_diff(ref: str, sub: str):
        if len(ref) != len(sub):
            return None

        return [(i, c2) for i, (c1, c2) in enumerate(zip(ref, sub)) if c1 != c2]

    @staticmethod
    def hamming_distance_from_diff(df1: str, df2: str):
        if len(df1) == 0:
            return len(df2)
        elif len(df2) == 0:
            return len(df1)

        i1, i2, n_same_position, n_same_pair = 0, 0, 0, 0
        while (i1 < len(df1) and i2 < len(df2)):
            if df1[i1][0] > df2[i2][0]:
                i2 += 1
            elif df1[i1][0] < df2[i2][0]:
                i1 += 1
            else:
                n_same_position += 1
                if df1[i1][1] == df2[i2][1]:
                    n_same_pair += 1
                i1 += 1
                i2 += 1

        return len(df1) + len(df2) - n_same_position - n_same_pair

    @staticmethod
    def levenshtein_distance(s1: str, s2: str) -> int:
        return editdistance.eval(s1, s2)

    @staticmethod
    def to_bitstream(dna_sequence: str):
        stream1 = 0x00
        stream2 = 0x00
        for nucl in dna_sequence:
            if nucl == 'A':
                pass
            elif nucl == 'G':
                stream1 += 1
            elif nucl == 'C':
                stream2 += 1
            elif nucl == 'T':
                stream1 += 1
                stream2 += 1
            else:
                return None
            stream1 = stream1 << 1
            stream2 = stream2 << 1

        return stream1, stream2

    @staticmethod
    def hamming_distance_from_bitstream(stream_pair1, stream_pair2):
        return bin((stream_pair1[0] ^ stream_pair2[0]) | (stream_pair1[1] ^ stream_pair2[1])).count('1')

    def translate(self, dna: str, codon_map: Dict[str, str] = None, **kwargs) -> str:
        if not codon_map:
            codon_map = self.standard_codon_map

        if codon_map == 'amber':
            codon_map = self.amber_codon_map

        aa = ''
        for i in range(0, len(dna), 3):
            triplet = dna[i:i + 3]
            try:
                aa += codon_map[triplet]
            except KeyError:
                if len(triplet) < 3:
                    aa += 'X'
                else:
                    aa += 'X'
            except Exception as e:
                print(e)

        return aa

    @staticmethod
    def RepresentsInt(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    def annotate_regions_SHN(self, chain_type: ChainType, sequence_list: List[str], annotation_list: List[str], **kwargs) -> TableData:
        """
        :param chain_type:
        :param sequence_list:
        :return:
        """

        # check sequence type
        # is_nucl = is_nucleotide_or_protein(sequence_list[0])
        #
        # if is_nucl == None or is_nucl == False:
        #     raise RuntimeError("parameter sequence_list is not a nucleotide sequence list")

        if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT,
                          ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_GAMMA,
                          ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT,
                          ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA,
                          ChainType.MOUSE_C57BL6_HEAVY]:
            tp_query = tempfile.NamedTemporaryFile('wt', suffix='.fasta', delete=False)
            if 'file_blast' in kwargs:
                file_blast = kwargs['file_blast']
            else:
                tp = tempfile.NamedTemporaryFile('wt')
                file_blast = tp.name
                tp.close()

            if 'domain_system' in kwargs:
                domain_system = kwargs['domain_system']
            else:
                domain_system = 'kabat'

            if 'blast_kwargs' in kwargs:
                blast_kwargs = kwargs['blast_kwargs']
            else:
                blast_kwargs = {}

            if 'seq_id_list' in kwargs:
                seq_id_list = kwargs['seq_id_list']
                for i, (seq, seq_id) in enumerate(zip(sequence_list, seq_id_list)):
                    tp_query.write(">%s\n" % seq_id)
                    tp_query.write(seq + '\n')

            else:
                for i, seq in enumerate(sequence_list):
                    tp_query.write(">seq_%d\n" % i)
                    tp_query.write(seq + '\n')

            tp_query.close()
            run_igblast(chain_type, tp_query.name, file_blast, domain_system, **blast_kwargs)
            print('[LOG] make %s done' % file_blast)

            result_list = self.parse_igblast(chain_type=chain_type,
                                             igblast_output=file_blast,
                                             domain_system=domain_system)

            # return result_list

            for sequence, blast_result in zip(sequence_list, result_list):
                aa_seq = self.translate(sequence)
                if len(sequence) % 3 == 0:
                    frame = 'In'
                else:
                    frame = 'Out'

                blast_result['query'] = sequence
                blast_result['frame'] = frame
                blast_result['stop_codon'] = '*' in aa_seq

                if 'FR4_from' in blast_result:
                    blast_result['FR4'] = sequence[blast_result['FR4_from']:blast_result['FR4_from'] + 33]

                for r in ['FR1', 'FR2', 'FR3', 'CDR1', 'CDR2']:
                    if ('%s_from' % r) in blast_result and ('%s_to' % r) in blast_result:
                        blast_result[r + '_NT'] = sequence[blast_result['%s_from' % r]:blast_result['%s_to' % r]]
                        blast_result[r + '_AA'] = self.translate(blast_result[r + '_NT'])

            newResultList = []
            for result in result_list:
                tempResultSet = {}
                for annotation in annotation_list:
                    try:
                        tempResultSet[annotation] = result[annotation]
                    except KeyError:
                        tempResultSet[annotation] = 'N/A'
                newResultList.append(tempResultSet)
            return newResultList
        else:
            raise RuntimeError('chain type error')

    def annotate_regions(self, chain_type: ChainType, sequence_list: List[str], annotation_list: List[str], **kwargs) -> TableData:
        """
        :param chain_type: 
        :param sequence_list: 
        :return: 
        """

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])

        if is_nucl == None or is_nucl == False:
            raise RuntimeError("parameter sequence_list is not a nucleotide sequence list")

        if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT, ChainType.CHICKEN_HEAVY,
                          ChainType.CHICKEN_LIGHT]:
            tp_query = tempfile.NamedTemporaryFile('wt', suffix='.fasta', delete=False)
            if 'file_blast' in kwargs:
                file_blast = kwargs['file_blast']
            else:
                tp = tempfile.NamedTemporaryFile('wt')
                file_blast = tp.name
                tp.close()

            if 'domain_system' in kwargs:
                domain_system = kwargs['domain_system']
            else:
                domain_system = 'kabat'

            if 'blast_kwargs' in kwargs:
                blast_kwargs = kwargs['blast_kwargs']
            else:
                blast_kwargs = {}

            if 'seq_id_list' in kwargs:
                seq_id_list = kwargs['seq_id_list']
                for i, (seq, seq_id) in enumerate(zip(sequence_list, seq_id_list)):
                    tp_query.write(">%s\n" % seq_id)
                    tp_query.write(seq + '\n')

            else:
                for i, seq in enumerate(sequence_list):
                    tp_query.write(">seq_%d\n" % i)
                    tp_query.write(seq + '\n')

            tp_query.close()
            run_igblast(chain_type, tp_query.name, file_blast, domain_system, **blast_kwargs)
            print('[LOG] make %s done' % file_blast)

            result_list = self.parse_igblast(chain_type=chain_type,
                                             igblast_output=file_blast,
                                             domain_system=domain_system)

            for sequence, blast_result in zip(sequence_list, result_list):
                aa_seq = self.translate(sequence)
                if len(sequence) % 3 == 0:
                    frame = 'In'
                else:
                    frame = 'Out'
                blast_result['query'] = sequence
                blast_result['frame'] = frame
                blast_result['stop codon'] = '*' in aa_seq
                if 'FR4 from' in blast_result:
                    blast_result['FR4'] = sequence[blast_result['FR4 from']:blast_result['FR4 from'] + 33]
                for r in ['FR1', 'FR2', 'FR3', 'CDR1', 'CDR2']:
                    if ('%s from' % r) in blast_result and ('%s to' % r) in blast_result:
                        blast_result[r] = sequence[blast_result['%s from' % r]:blast_result['%s to' % r]]

            new_result_list = []
            for result in result_list:
                temp_new_result = []
                for annotation in annotation_list:
                    try:
                      temp_new_result.append(result[annotation])
                    except KeyError:
                        temp_new_result.append('N/A')
                new_result_list.append(temp_new_result)
            return new_result_list
        else:
            raise RuntimeError('chain type error')

    def annotate_regions_chicken(self, seq_list, v_type, result_dir='/home/IP-team/data/blast_db/blast_temp/'):
        query_dir = os.path.join(result_dir, 'query')
        out_dir = os.path.join(result_dir, 'out')

        if os.path.exists(query_dir) is not True:
            os.makedirs(query_dir)

        if os.path.exists(out_dir) is not True:
            os.makedirs(out_dir)

        v_type = v_type.lower()

        with open(os.path.join(query_dir, 'query.fasta'), 'w') as fasta_writer:
            for i, e_seq in enumerate(seq_list):
                fasta_writer.write('>'+str(i)+'\n')
                fasta_writer.write(e_seq+'\n')

        db_dir = '/home/IP-team/data/blast_db/fr_db_chicken_white_leghorn'

        region_list = ['fr1','fr2','fr3','fr4']

        # Variables used in region extraction & Initialization
        ext_index = []
        cdr_list = []

        for i in range(len(seq_list)):
            temp = {'fr1_e': None, 'fr2_s': None, 'fr2_e': None, 'fr3_s': None, 'fr3_e': None, 'fr4_s': None}
            ext_index.append(temp)

        # Parameters for blast_run
        E_VAL = 5
        NUM_ALIGN = 10
        THREADS = 20
        FORMAT_EXE = '/Tools/ncbi-blast-2.7.1+/bin/blastn'

        for region in region_list:
            db_name = '_'.join([v_type, region])
            db_file = os.path.join(db_dir, db_name)
            out_name = 'query_' + db_name + '_blasted.txt'
            out_file = os.path.join(out_dir, out_name)

            print('Blast run for ' + db_name + ' is started.')

            # Note that word_size = 5
            cline = NcbiblastnCommandline(num_threads=THREADS, query=os.path.join(query_dir, 'query.fasta'), db=db_file, evalue=0.1** E_VAL, \
                                          out=out_file, gapopen=1, gapextend=2, word_size=5, num_descriptions=NUM_ALIGN, num_alignments=NUM_ALIGN)

            os.system(FORMAT_EXE + str(cline)[len('blastn'):])

            print('Blast run for ' + db_name + ' is completed.\n')

            handle = open(out_file, 'r')
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            for i, each in enumerate(iterator):
                if len(each.alignments) == 0:
                    continue

                for alignment in each.alignments:
                    hsps = alignment.hsps[0]

                    s_point = hsps.query_start
                    e_point = hsps.query_end

                    if region == 'fr1':
                        if ext_index[i]['fr1_e'] == None and hsps.sbjct_end == alignment.length:
                            ext_index[i]['fr1_e'] = e_point
                            break
                    elif region == 'fr2':
                        if ext_index[i]['fr2_s'] == None and hsps.sbjct_start == 1:
                            ext_index[i]['fr2_s'] = s_point - 1
                        if ext_index[i]['fr2_e'] == None and hsps.sbjct_end == alignment.length:
                            ext_index[i]['fr2_e'] = e_point
                        if ext_index[i]['fr2_s'] != None and ext_index[i]['fr2_e'] != None:
                            break
                    elif region == 'fr3':
                        if ext_index[i]['fr3_s'] == None and hsps.sbjct_start == 1:
                            ext_index[i]['fr3_s'] = s_point - 1
                        if ext_index[i]['fr3_e'] == None and hsps.sbjct_end == alignment.length:
                            ext_index[i]['fr3_e'] = e_point
                        if ext_index[i]['fr3_s'] != None and ext_index[i]['fr3_e'] != None:
                            break
                    elif region == 'fr4':
                        if ext_index[i]['fr4_s'] == None and hsps.sbjct_start == 1:
                            ext_index[i]['fr4_s'] = s_point - 1
                            break

        for seq, ext in zip(seq_list, ext_index):
            if ext['fr1_e'] == None or ext['fr2_s'] == None:
                cdr1 = 'E/F'
            else:
                cdr1 = seq[ext['fr1_e']:ext['fr2_s']]

            if ext['fr2_e'] == None or ext['fr3_s'] == None:
                cdr2 = 'E/F'
            else:
                cdr2 = seq[ext['fr2_e']:ext['fr3_s']]

            if ext['fr3_e'] == None or ext['fr4_s'] == None:
                cdr3 = 'E/F'
            else:
                cdr3 = seq[ext['fr3_e']:ext['fr4_s']]

            cdr_list.append({'cdr1':cdr1, 'cdr2':cdr2, 'cdr3':cdr3})

        fail_list = []
        for i, e_cdr in enumerate(cdr_list):
            if e_cdr['cdr3'] == 'E/F':
                fail_list.append(i)

        return cdr_list, fail_list

    def blast_run(self, query_list: list(), db_list: list(), m_type: str, blast_dir: str = '/home/IP-team/data/blast_temp/', **kwargs):
        if 'file_blast' in kwargs:
            file_blast = kwargs['file_blast']

        if 'num_threads' in kwargs:
            num_threads = kwargs['num_threads']
        else:
            num_threads = 10

        if 'evalue' in kwargs:
            evalue = kwargs['evalue']
        else:
            evalue = 0.1**5

        if 'gapopen' in kwargs:
            gapopen = kwargs['gapopen']
        else:
            gapopen = 1

        if 'gapextend' in kwargs:
            gapextend = kwargs['gapextend']
        else:
            gapextend = 2

        if 'word_size' in kwargs:
            word_size = kwargs['word_size']
        else:
            word_size = 10

        if 'num_alignments' in kwargs:
            num_alignments = kwargs['num_alignments']
        else:
            num_alignments = 1

        if 'result_file_name' in kwargs:
            result_file_name = kwargs['result_file_name']
        else:
            result_file_name = 'blasted.txt'

        query_dir = os.path.join(blast_dir, 'query')
        db_dir = os.path.join(blast_dir, 'db')
        result_dir = os.path.join(blast_dir, 'result')

        if not os.path.exists(query_dir):
            os.makedirs(query_dir)

        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        result_file = os.path.join(result_dir, result_file_name)

        with open(os.path.join(query_dir, 'query.fasta'), 'w') as fasta_writer:
            for i, e_seq in enumerate(query_list):
                fasta_writer.write('>'+str(i)+'\n')
                fasta_writer.write(e_seq+'\n')

        with open(os.path.join(db_dir, 'db.fasta'), 'w') as fasta_writer:
            for i, e_seq in enumerate(db_list):
                fasta_writer.write('>'+str(i)+'\n')
                fasta_writer.write(e_seq+'\n')

        if m_type.upper() == 'AA':
            # construct db
            format_cmd = '%s -in %s -dbtype prot -input_type fasta -out %s' % ('/Tools/ncbi-blast-2.7.1+/bin/makeblastdb', os.path.join(db_dir, 'db.fasta'), os.path.join(db_dir, 'db'))
            os.system(format_cmd)

            cline = NcbiblastpCommandline(num_threads= num_threads, query= os.path.join(query_dir, 'query.fasta'), db= os.path.join(db_dir, 'db'), out= result_file,
                                          num_descriptions=num_alignments, num_alignments= num_alignments, evalue=evalue)

            format_exe = '/Tools/ncbi-blast-2.7.1+/bin/blastp'

            os.system(format_exe + str(cline)[len('blastp'):])

            print('Blastp completed.\n')

            handle = open(result_file, 'r')
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            resultList = []

            for i, each in enumerate(iterator):
                numAlignments = 0
                alignmentList = []

                if len(each.alignments) == 0:
                    tempSet = {'num_alignments':numAlignments, 'alignments':alignmentList}
                    resultList.append(tempSet)
                else:
                    for alignment in each.alignments:
                        sbjctIndex = int(alignment.title[1:].strip())
                        hsps = alignment.hsps[0]
                        numAlignments += 1
                        queryStart = hsps.query_start
                        queryEnd = hsps.query_end
                        sbjctStart = hsps.sbjct_start
                        sbjctEnd = hsps.sbjct_end

                        matched = hsps.identities[0]
                        whole = hsps.identities[1]

                        numMismatches = hsps.identities[1] - hsps.identities[0]
                        mismatchList = []
                        for i, eMatch in enumerate(hsps.match):
                            if eMatch == ' ':
                                queryLoc = i + queryStart
                                sbjctLoc = i + sbjctStart
                                queryMismatch = hsps.query[i]
                                sbjctMismatch = hsps.sbjct[i]
                                mismatchList.append({'query':queryMismatch, 'sbjct': sbjctMismatch, 'query_loc': queryLoc, 'sbjct_loc': sbjctLoc})
                        alignmentList.append({'matched':matched, 'whole':whole, 'sbjct_index': sbjctIndex, 'query_start':queryStart, 'query_end': queryEnd, 'sbjct_start': sbjctStart, 'sbjct_end': sbjctEnd, 'num_mismatches': numMismatches, 'mismatches': mismatchList})
                    tempSet = {'num_alignments':numAlignments, 'alignments': alignmentList}
                    resultList.append(tempSet)
            return resultList

        if m_type.upper() == 'NT':
            # construct db
            format_cmd = '%s -in %s -dbtype nucl -input_type fasta -out %s' % ('/Tools/ncbi-blast-2.7.1+/bin/makeblastdb', os.path.join(db_dir, 'db.fasta'), os.path.join(db_dir, 'db'))
            os.system(format_cmd)

            cline = NcbiblastnCommandline(num_threads=num_threads, query=os.path.join(query_dir, 'query.fasta'), db=os.path.join(db_dir, 'db'),
                                          evalue=evalue, out=result_file, gapopen=gapopen, gapextend=gapextend, word_size=word_size,
                                          num_descriptions=num_alignments, num_alignments=num_alignments)


            format_exe = '/Tools/ncbi-blast-2.7.1+/bin/blastn'

            os.system(format_exe + str(cline)[len('blastn'):])

            print('Blastn completed.\n')

            handle = open(result_file, 'r')
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            resultList = []

            for i, each in enumerate(iterator):
                numAlignments = 0
                alignmentList = []

                if len(each.alignments) == 0:
                    tempSet = {'num_alignments': numAlignments, 'alignments': alignmentList}
                    resultList.append(tempSet)
                else:
                    for alignment in each.alignments:
                        sbjctIndex = int(alignment.title[1:].strip())
                        hsps = alignment.hsps[0]
                        numAlignments += 1
                        queryStart = hsps.query_start
                        queryEnd = hsps.query_end
                        sbjctStart = hsps.sbjct_start
                        sbjctEnd = hsps.sbjct_end

                        matched = hsps.identities[0]
                        whole = hsps.identities[1]

                        numMismatches = hsps.identities[1] - hsps.identities[0]
                        mismatchList = []
                        for i, eMatch in enumerate(hsps.match):
                            if eMatch == ' ':
                                queryLoc = i + queryStart
                                sbjctLoc = i + sbjctStart
                                queryMismatch = hsps.query[i]
                                sbjctMismatch = hsps.sbjct[i]
                                mismatchList.append(
                                    {'query': queryMismatch, 'sbjct': sbjctMismatch, 'query_loc': queryLoc,
                                     'sbjct_loc': sbjctLoc})
                        alignmentList.append(
                            {'matched': matched, 'whole': whole, 'sbjct_index': sbjctIndex, 'query_start': queryStart,
                             'query_end': queryEnd, 'sbjct_start': sbjctStart, 'sbjct_end': sbjctEnd,
                             'num_mismatches': numMismatches, 'mismatches': mismatchList})
                    tempSet = {'num_alignments': numAlignments, 'alignments': alignmentList}
                    resultList.append(tempSet)
            return resultList

    def annotate_regions_chicken_loose(self, seq_list, v_type, result_dir='/home/IP-team/data/blast_db/blast_temp/'):
        query_dir = os.path.join(result_dir, 'query')
        out_dir = os.path.join(result_dir, 'out')

        if os.path.exists(query_dir) is not True:
            os.makedirs(query_dir)

        if os.path.exists(out_dir) is not True:
            os.makedirs(out_dir)

        v_type = v_type.lower()

        with open(os.path.join(query_dir, 'query.fasta'), 'w') as fasta_writer:
            for i, e_seq in enumerate(seq_list):
                fasta_writer.write('>'+str(i)+'\n')
                fasta_writer.write(e_seq+'\n')

        db_dir = '/home/IP-team/data/blast_db/fr_db_chicken_white_leghorn'

        region_list = ['fr1','fr2','fr3','fr4']

        # Variables used in region extraction & Initialization
        ext_index = []
        anno_list = []

        for i in range(len(seq_list)):
            temp = {'fr1_s':None,'fr1_e': None, 'fr2_s': None, 'fr2_e': None, 'fr3_s': None, 'fr3_e': None, 'fr4_s': None,'fr4_e':None}
            ext_index.append(temp)

        # Parameters for blast_run
        E_VAL = 5
        NUM_ALIGN = 1
        THREADS = 20
        FORMAT_EXE = '/Tools/ncbi-blast-2.7.1+/bin/blastn'

        for region in region_list:
            db_name = '_'.join([v_type, region])
            db_file = os.path.join(db_dir, db_name)
            out_name = 'query_' + db_name + '_blasted.txt'
            out_file = os.path.join(out_dir, out_name)

            print('Blast run for ' + db_name + ' is started.')

            # Note that word_size = 10
            cline = NcbiblastnCommandline(num_threads=THREADS, query=os.path.join(query_dir, 'query.fasta'), db=db_file, evalue=0.1** E_VAL, \
                                          out=out_file, gapopen=1, gapextend=2, word_size=10, num_descriptions=NUM_ALIGN, num_alignments=NUM_ALIGN)

            os.system(FORMAT_EXE + str(cline)[len('blastn'):])

            print('Blast run for ' + db_name + ' is completed.\n')

            handle = open(out_file, 'r')
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            for i, each in enumerate(iterator):
                if len(each.alignments) == 0:
                    continue

                for alignment in each.alignments:
                    hsps = alignment.hsps[0]

                    s_point = hsps.query_start - hsps.sbjct_start
                    e_point = hsps.query_end + alignment.length  - hsps.sbjct_end

                    if region == 'fr1':
                        ext_index[i]['fr1_s'] = s_point
                        ext_index[i]['fr1_e'] = e_point
                        break
                    elif region == 'fr2':
                        ext_index[i]['fr2_s'] = s_point
                        ext_index[i]['fr2_e'] = e_point
                        break
                    elif region == 'fr3':
                        ext_index[i]['fr3_s'] = s_point
                        ext_index[i]['fr3_e'] = e_point
                        break
                    elif region == 'fr4':
                        ext_index[i]['fr4_s'] = s_point
                        ext_index[i]['fr4_e'] = e_point
                        break

        for seq, ext in zip(seq_list, ext_index):
            if ext['fr1_e'] == None or ext['fr2_s'] == None:
                cdr1 = 'E/F'
            else:
                cdr1 = seq[ext['fr1_e']:ext['fr2_s']]

            if ext['fr2_e'] == None or ext['fr3_s'] == None:
                cdr2 = 'E/F'
            else:
                cdr2 = seq[ext['fr2_e']:ext['fr3_s']]

            if ext['fr3_e'] == None or ext['fr4_s'] == None:
                cdr3 = 'E/F'
            else:
                cdr3 = seq[ext['fr3_e']:ext['fr4_s']]

            if ext['fr1_s'] == None or ext['fr4_e'] == None:
                full = 'E/F'
            else:
                full = seq[ext['fr1_s']:ext['fr4_e']]

            if ext['fr1_e'] == None or ext['fr4_s'] == None:
                fullWoFr = 'E/F'
            else:
                fullWoFr = seq[ext['fr1_e']:ext['fr4_s']]

            anno_list.append({'full':full, 'full_wo_fr':fullWoFr, 'cdr1':cdr1, 'cdr2':cdr2, 'cdr3':cdr3})

        return anno_list

    def parse_igblast(self, chain_type: ChainType, igblast_output: str, **kwargs):
        result_list = []
        result = None

        if 'domain_system' in kwargs:
            domain_system = kwargs['domain_system']
        else:
            domain_system = 'kabat'

        with open(igblast_output, 'r') as igb_handle:
            if chain_type in [ChainType.HUMAN_LIGHT, ChainType.HUMAN_HEAVY]:
                while True:
                    line = igb_handle.readline()
                    if not line:
                        break

                    if line[:len('Query= ')] == 'Query= ':
                        if result is not None:
                            result_list.append(result)
                        result = {'query': line[len('Query= '):-1], 'hit': False, 'V gene': '', 'J gene': '',
                                  'CDR3': '-'}
                        if chain_type == ChainType.HUMAN_HEAVY:
                            result['D gene'] = ''
                    elif 'No hits found' in line:
                        result['hit'] = False
                    elif line[:len('Sub-region sequence details')] == 'Sub-region sequence details':
                        next_line_sp = igb_handle.readline().split('\t')
                        if next_line_sp[0] == 'CDR3':
                            result['hit'] = True
                            result['CDR3'] = next_line_sp[1]
                            result['FR4 from'] = int(next_line_sp[4])
                    elif line[:len('Alignment summary')] == 'Alignment summary':
                        while True:
                            next_line_sp = igb_handle.readline().split('\t')
                            if domain_system == 'imgt':
                                region = next_line_sp[0].replace('-IMGT', '')
                            else:
                                region = next_line_sp[0]
                            if region == 'Total':
                                break
                            if region == 'FR1':
                                fr1_start = int(next_line_sp[1]) - 1
                                fr1_end = int(next_line_sp[2])
                                if self.RepresentsInt(next_line_sp[6]):
                                    fr1_gap = int(next_line_sp[6])
                                else:
                                    fr1_gap = 0
                                while ((fr1_end - fr1_start + fr1_gap) % 3 != 0):
                                    fr1_start += 1
                                result['FR1 from'] = fr1_start
                                result['FR1 to'] = fr1_end
                            else:
                                result[region + ' from'] = int(next_line_sp[1]) - 1
                                result[region + ' to'] = int(next_line_sp[2])
                    elif line[:len('V-(D)-J rearrangement summary')] == 'V-(D)-J rearrangement summary':
                        next_line_sp = igb_handle.readline().split('\t')
                        if chain_type == ChainType.HUMAN_HEAVY:
                            result['V gene'] = next_line_sp[0]
                            result['D gene'] = next_line_sp[1]
                            result['J gene'] = next_line_sp[2]
                            result['V-J frame'] = next_line_sp[5].split('-')[0]
                        elif chain_type == ChainType.HUMAN_LIGHT:
                            result['V gene'] = next_line_sp[0]
                            result['J gene'] = next_line_sp[1]
                            result['V-J frame'] = next_line_sp[4].split('-')[0]
            elif chain_type in [ChainType.CHICKEN_LIGHT, ChainType.CHICKEN_HEAVY]:
                while True:
                    line = igb_handle.readline()
                    if not line:
                        break

                    if line[:len('Query= ')] == 'Query= ':
                        if result is not None:
                            result_list.append(result)
                        result = {'query': line[len('Query= '):-1], 'hit': False, 'CDR3': ''}
                    elif 'No hits found' in line:
                        result['hit'] = False
                    elif line[:len('V-(D)-J junction details')] == 'V-(D)-J junction details':
                        next_line_sp = igb_handle.readline().split('\t')
                        if chain_type == ChainType.CHICKEN_HEAVY:
                            result['CDR3'] = next_line_sp[3]
                        else:
                            result['CDR3'] = next_line_sp[1]
                        if result['CDR3'] != '':
                            result['hit'] = True

        if result is not None:
            result_list.append(result)

        return result_list

    def numbering_CDR3(self, chain_type: ChainType, cdr3_list: List[str], alphabet_with_number: bool = True) -> Tuple[
        TableHeader, TableData]:

        # check sequence type
        is_nucl = is_nucleotide_or_protein(cdr3_list[0])
        if is_nucl == None:
            raise RuntimeError("parameter cdr3_list unknown")

        if chain_type in (ChainType.HUMAN_LIGHT, ChainType.CHICKEN_LIGHT):
            CDR3_numbering = ['24', '25', '26', '27', '28', '29', '30', 'A', '31', '32', '33', '34']
            len_end = 4
            cdr3_key = 'LCDR3'
        elif chain_type in (ChainType.HUMAN_HEAVY, ChainType.CHICKEN_HEAVY):
            CDR3_numbering = ['95', '96', '97', '98', '99', '100', 'A', '101', '102']
            len_end = 2
            cdr3_key = 'HCDR3'
        else:
            raise RuntimeError("parameter chain_type error")
        i_A = CDR3_numbering.index('A')
        if alphabet_with_number:
            last_num = CDR3_numbering[i_A - 1]
        else:
            last_num = ''
        CDR3_numbering[i_A] = last_num + CDR3_numbering[i_A]

        if is_nucl:
            max_len = int(max([len(s) / 3 for s in cdr3_list]))
        else:
            max_len = int(max([len(s) for s in cdr3_list]))
        init_len = len(CDR3_numbering)
        if max_len > init_len:
            CDR3_numbering = CDR3_numbering[:i_A + 1] + [last_num + chr(ord("A") + i + 1) for i in
                                                         range(max_len - init_len)] + CDR3_numbering[i_A + 1:]

        return_data = []
        for cdr3_str in cdr3_list:
            if cdr3_str == '-' or cdr3_str == 'X' or cdr3_str == '':
                d = {cdr3_key: '-', 'query': cdr3_str}
                for n in CDR3_numbering:
                    d[n] = ''
                return_data.append(d)
                continue

            if is_nucl:
                cdr3_str_aa = self.translate(cdr3_str)
            else:
                cdr3_str_aa = cdr3_str
            cdr3_len = len(cdr3_str_aa)
            d = {cdr3_key: cdr3_str_aa, 'query': cdr3_str_aa}
            for n in CDR3_numbering:
                d[n] = ''

            if cdr3_len == 1:
                d[CDR3_numbering[0]] = cdr3_str_aa[0]
            else:
                d[CDR3_numbering[-1]] = cdr3_str_aa[-1]
                len_end_this = min(len_end, cdr3_len - 1)
                for i in range(-1, -1 * len_end_this - 1, -1):
                    d[CDR3_numbering[i]] = cdr3_str_aa[i]
                for i in range(cdr3_len - len_end_this):
                    d[CDR3_numbering[i]] = cdr3_str_aa[i]
            return_data.append(d)

        return CDR3_numbering, return_data

    def numbering_all(self, chain_type: ChainType, data: TableData,
                      format_empty_residue='{bg_color:gray}') -> Tuple[TableHeader, List[TableHeader], List[TableData]]:

        if chain_type in (ChainType.HUMAN_LIGHT, ChainType.CHICKEN_LIGHT):
            numbering = {
                'FR1': ['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9',
                        'L10', 'L11', 'L12', 'L13', 'L14', 'L15', 'L16', 'L17', 'L18', 'L19',
                        'L20', 'L21', 'L22', 'L23', ],
                'CDR1': ['L24', 'L25', 'L26', 'L27', 'L27A', 'L28', 'L29', 'L30', 'L31', 'L32', 'L33',
                         'L34'],
                'FR2': ['L35', 'L36', 'L37', 'L38', 'L39', 'L40', 'L41', 'L42', 'L43', 'L44',
                        'L45', 'L46', 'L47', 'L48', 'L49'],
                'CDR2': ['L50', 'L51', 'L52', 'L53', 'L54', 'L55', 'L56'],
                'FR3': ['L57', 'L58', 'L59', 'L60', 'L61', 'L62', 'L63', 'L64', 'L65', 'L66',
                        'L67', 'L68', 'L69', 'L70', 'L71', 'L72', 'L73', 'L74', 'L75', 'L76',
                        'L77', 'L78', 'L79', 'L80', 'L81', 'L82', 'L83', 'L84', 'L85', 'L86',
                        'L87', 'L88', ],
                'CDR3': ['L89', 'L90', 'L91', 'L92', 'L93', 'L94', 'L95', 'L95A', 'L96', 'L97'],
                'FR4': ['L98', 'L99', 'L100', 'L101', 'L102', 'L103', 'L104', 'L105', 'L106', 'L106A', 'L107',
                        'L108', 'L109']
            }
        elif chain_type in (ChainType.HUMAN_HEAVY, ChainType.CHICKEN_HEAVY):
            numbering = {
                'FR1': ['H0', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9',
                        'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'H17', 'H18', 'H19',
                        'H20', 'H21', 'H22', 'H23', 'H24', 'H25', 'H26', 'H27', 'H28', 'H29',
                        'H30'],
                'CDR1': ['H31', 'H32', 'H33', 'H34', 'H35', 'H35A'],
                'FR2': ['H36', 'H37', 'H38', 'H39', 'H40', 'H41', 'H42', 'H43', 'H44', 'H45',
                        'H46', 'H47', 'H48', 'H49'],
                'CDR2': ['H50', 'H51', 'H52', 'H52A', 'H53', 'H54', 'H55', 'H56', 'H57', 'H58', 'H59',
                         'H60', 'H61', 'H62', 'H63', 'H64', 'H65'],
                'FR3': ['H66', 'H67', 'H68', 'H69', 'H70', 'H71', 'H72', 'H73', 'H74', 'H75',
                        'H76', 'H77', 'H78', 'H79', 'H80', 'H81', 'H82', 'H82A', 'H83', 'H84', 'H85',
                        'H86', 'H87', 'H88', 'H89', 'H90', 'H91', 'H92', 'H93', 'H94'],
                'CDR3': ['H95', 'H96', 'H97', 'H98', 'H99', 'H100', 'H100A', 'H101', 'H102'],
                'FR4': ['H103', 'H104', 'H105', 'H106', 'H107', 'H108', 'H109', 'H110', 'H111', 'H112',
                        'H113']
            }
        else:
            raise RuntimeError("parameter chain_type error")

        regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
        list_numbering_data = []
        list_numbering_header = [numbering[r] for r in regions]
        for i_r, (r, numbernig_header) in enumerate(zip(regions, list_numbering_header)):
            numbering_data = []

            r_max_length = 0
            for d in data:
                d_num = {}
                if r not in d:
                    numbering_data.append(d_num)
                    continue

                is_nucl = is_nucleotide_or_protein(d[r])
                if is_nucl == None:
                    raise RuntimeError("unknown sequence type")
                if is_nucl:
                    r_seq_aa = self.translate(d[r])
                else:
                    r_seq_aa = d[r]
                r_max_length = max(r_max_length, len(r_seq_aa))

                last_index_alpha = -1
                for i in range(len(numbernig_header)):
                    if numbernig_header[-i - 1][-1].isalpha():
                        last_index_alpha = len(numbernig_header) - i - 1
                        break
                # if length is variable
                if r_max_length > len(numbernig_header):
                    if last_index_alpha != -1:
                        deficient_length = r_max_length - len(numbernig_header)
                        last_alpha = numbernig_header[last_index_alpha][-1]
                        numbernig_header = numbernig_header[:last_index_alpha + 1] + [
                            numbernig_header[last_index_alpha][:-1] + chr(ord(last_alpha) + j + 1) for j in
                            range(deficient_length)] + numbernig_header[last_index_alpha + 1:]
                        last_index_alpha += deficient_length
                        list_numbering_header[i_r] = numbernig_header

                if r == 'FR1':
                    stick_left = 0
                elif r == 'FR4':
                    stick_left = len(r_seq_aa)
                else:
                    # if do not have alphabet part
                    if last_index_alpha == -1:
                        stick_left = len(r_seq_aa)
                    else:
                        stick_left = len(r_seq_aa) - (len(numbernig_header) - last_index_alpha - 1)

                    if stick_left < 0:
                        stick_left = 0

                for i, n in enumerate(numbernig_header):
                    if i == stick_left:
                        break
                    d_num[n] = r_seq_aa[i]
                for j, n in enumerate(numbernig_header[::-1]):
                    if len(r_seq_aa) - stick_left == j:
                        break
                    d_num[n] = r_seq_aa[-j - 1]

                numbering_data.append(d_num)

            list_numbering_data.append(numbering_data)

        # fill the blank
        for r, numbernig_header, numbering_data in zip(regions, list_numbering_header, list_numbering_data):
            for d_num in numbering_data:
                for n in numbernig_header:
                    if n not in d_num:
                        d_num[n] = format_empty_residue + '-'

        return regions, list_numbering_header, list_numbering_data

    def add_ptm_count_and_color(self, aa_seq_list=List[str], prefix_residue: str = 'R', **kwargs) -> Tuple[
        TableHeader, TableHeader, TableData]:
        """
        additional_ptm = [('D', '(D)'), ('E', '(E)'), ('D+E', '([DE])'), ('R', '(R)'), ('K', '(K)'), ('R+K', '([RK])')]

        :param header: 
        :param data: 
        :param col_aa_seq: 
        :param col_residue: Typically, ('24', '34') for VL and ('95', '102') for VH.
                            If this value is None, just add count info and do not mark color.
        :return: 
        """

        ret_data = []
        if 'aa_residues' in kwargs and 'resi_header' in kwargs:
            ret_data = kwargs['aa_residues']
            resi_header = kwargs['resi_header']
        else:
            lengths = [len(seq) for seq in aa_seq_list]
            max_length = int(max(lengths))
            for seq in aa_seq_list:
                d = {}
                d['query'] = seq
                for n in range(max_length):
                    if n < len(seq):
                        d[prefix_residue + str(n + 1)] = seq[n]
                    else:
                        d[prefix_residue + str(n + 1)] = ''
                ret_data.append(d)
            resi_header = [prefix_residue + str(n + 1) for n in range(max_length)]

        def mark_ptm(_d, _col_letters, _ptm_name, mark_str, match_list):
            _d[_ptm_name] = len(match_list)
            if len(_col_letters) > 0:
                for m in match_list:
                    for i in range(*m.span()):
                        _d[_col_letters[i]] = mark_str + _d[_col_letters[i]]

        default_ptms = ['Deamidation', 'Aspartate_Isomerization', 'N Glycosylation', 'Cleavage', 'Oxidation',
                        'Free Cysteine', 'Sulfation', 'Methylation']

        finder = PTMFinder()
        for d in ret_data:
            seq = d['query']
            col_letters = []
            for h in resi_header:
                try:
                    if d[h] != '' and detach_format(d[h]) == seq[len(col_letters)]:
                        col_letters.append(h)
                except Exception as ex:
                    print(ex)

            mark_ptm(d, col_letters, default_ptms[0], "{bg_color:blue}", finder.find_diamidation(seq))
            mark_ptm(d, col_letters, default_ptms[1], "{bg_color:cyan}", finder.find_aspartate_isomerization(seq))
            mark_ptm(d, col_letters, default_ptms[2], "{bg_color:gray}", finder.find_n_glycosylation(seq))
            mark_ptm(d, col_letters, default_ptms[3], "{bg_color:green}", finder.find_cleavage(seq))
            mark_ptm(d, col_letters, default_ptms[4], "{bg_color:lime}", finder.find_oxidation(seq))
            mark_ptm(d, col_letters, default_ptms[5], "{bg_color:magenta}", finder.find_free_cysteine(seq))
            mark_ptm(d, col_letters, default_ptms[6], "{bg_color:orange}", finder.find_sulfation(seq))
            mark_ptm(d, col_letters, default_ptms[7], "{bg_color:yellow}", finder.find_methylation(seq))

        return default_ptms, resi_header, ret_data

    def get_ptm_count(self, aa_seq_list=List[str]) -> Tuple[TableHeader, TableData]:
        default_ptms = ['Deamidation', 'Aspartate_Isomerization', 'N Glycosylation', 'Cleavage', 'Oxidation',
                        'Free Cysteine', 'Sulfation', 'Methylation']

        ret_data = []
        finder = PTMFinder()
        for seq in aa_seq_list:
            d = {'query': seq,
                 default_ptms[0]: len(finder.find_diamidation(seq)),
                 default_ptms[1]: len(finder.find_aspartate_isomerization(seq)),
                 default_ptms[2]: len(finder.find_n_glycosylation(seq)),
                 default_ptms[3]: len(finder.find_cleavage(seq)),
                 default_ptms[4]: len(finder.find_oxidation(seq)),
                 default_ptms[5]: len(finder.find_free_cysteine(seq)),
                 default_ptms[6]: len(finder.find_sulfation(seq)),
                 default_ptms[7]: len(finder.find_methylation(seq)),
                 }
            ret_data.append(d)

        return default_ptms, ret_data

    def get_pi_value(self, aa_seq: str) -> float:
        putil = ProteinAnalysis(aa_seq)
        return putil.isoelectric_point()

    def get_gravy_value(self, aa_seq: str) -> float:
        putil = ProteinAnalysis(aa_seq)
        return putil.gravy()


class PhyloWorker(Process):
    def __init__(self, in_queue, out_queue, run_type):
        Process.__init__(self)

        self.in_queue = in_queue
        self.out_queue = out_queue
        self.run_type = run_type

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            target = self.in_queue.get()
            if target is None:
                self.in_queue.task_done()
                break
            if self.run_type == 'extract_edges_phylotree':
                self.extract_edges_phylotree(**target)
            elif self.run_type == 'build_tree':
                self.build_tree(**target)
            self.in_queue.task_done()

        logging.info("Finished run")

    def extract_edges_phylotree(self, col_group_id: str, group_name: str, msa_file: str):
        tree_output = msa_file.replace('msa', 'tree').replace('.fa', '.nex')
        is_exist = False
        if Path(tree_output).is_file() and os.stat(str(tree_output)).st_size > 0:
            tree = Phylo.read(tree_output, 'nexus')
            is_exist = True
        else:
            align = AlignIO.read(msa_file, "fasta")
            calculator = DistanceCalculator('blastn')
            dm = calculator.get_distance(align)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)

        edges = []

        def get_clade_name(clade):
            return group_name + '_' + clade.name if 'Inner' in clade.name else clade.name

        def get_edge_from_clade(clade, is_mother):
            if len(clade.clades) == 0:
                return
            edges.append({
                'from': get_clade_name(clade),
                'to': get_clade_name(clade.clades[0]),
                'length': clade.clades[0].branch_length,
                'root': is_mother,
                col_group_id: group_name,
            })
            edges.append({
                'from': get_clade_name(clade),
                'to': get_clade_name(clade.clades[1]),
                'length': clade.clades[1].branch_length,
                'root': 0,
                col_group_id: group_name,
            })
            get_edge_from_clade(clade.clades[0], 0)
            get_edge_from_clade(clade.clades[1], 0)

        get_edge_from_clade(tree.clade, 1)

        if not is_exist:
            Phylo.write(tree, msa_file.replace('msa', 'tree').replace('.fa', '.nex'), 'nexus')

        self.out_queue.put([group_name, tree.clade.name, edges])

    def build_tree(self, aln_file: str, distance_model: str, **kwargs):
        align = AlignIO.read(aln_file, "fasta")
        calculator = DistanceCalculator(distance_model)
        dm = calculator.get_distance(align)
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)
        self.out_queue.put({'tree': tree, **kwargs})

    @staticmethod
    def cut_tree(tree: Tree, branch_length_thres: float) -> List[Clade]:
        subtree_list = []

        def _dfs_cut_tree(_clade: Clade, _subtree_list: List[Clade]):
            new_children = []
            for i, _child in enumerate(_clade.clades):
                if _child.branch_length > branch_length_thres:
                    _child.branch_length = 0
                    _subtree_list.append(_child)
                else:
                    new_children.append(_child)
                _dfs_cut_tree(_child, _subtree_list)
            _clade.clades = new_children

        subtree_list.append(tree.clade)
        _dfs_cut_tree(tree.clade, subtree_list)

        return subtree_list

    @staticmethod
    def sort_sequence_by_phylo(sequence_batch, tree, out_file):
        map_reads = {}
        for rec in sequence_batch:
            map_reads[rec.id] = rec

        def dps_tree(clade, out_list):
            if len(clade.clades) == 0:
                out_list.append(clade.name)
            else:
                dps_tree(clade.clades[0], out_list)
                dps_tree(clade.clades[1], out_list)
            return

        sorted_id_list = []
        dps_tree(tree.clade, sorted_id_list)

        records_to_out = []
        for id in sorted_id_list:
            records_to_out.append(map_reads[id])

        SeqIO.write(records_to_out, out_file, 'fasta')


def insert_into_pComb3x(full_vl_seq: str, full_vh_seq: str, linker_seq: str) -> str:
    vector_seq_file = r'C:\Users\hhjunny\Lib\Pcomb3X_Vector.fa'
    before_insert = 'CCAGGCGGCC'
    after_insert = 'ACTAGTGGCC'
    vector = SeqIO.read(vector_seq_file, 'fasta')

    before_insert_full = str(vector.seq)[:str(vector.seq).index(before_insert) + len(before_insert)]
    after_insert_full = str(vector.seq)[str(vector.seq).index(after_insert):]
    return before_insert_full + full_vl_seq + linker_seq + full_vh_seq + after_insert_full
