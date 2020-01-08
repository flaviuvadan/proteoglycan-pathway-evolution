import os

import numpy as np
from rpy2 import robjects
from rpy2.robjects import packages, vectors


class Loader:
    """ Responsible for loading in organisms data for constructing dnds visualizations """

    GENE_NAME_IDX = 0

    ORG_CLASS_IDX = 0
    ORG_NAME_IDX = 1

    UNACCEPTED_CHARS = ["U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]

    def __init__(self):
        self.alignments = self._parse_organisms_alignments()
        self.significant_orgs = self._get_significant_orgs()
        self.get_dnds = self._setup_R_access()
        self.scores = self._get_dnds_scores()

    def _parse_organisms_alignments(self):
        """ Reads in the organisms alignments """
        dir_path = os.path.join(os.getcwd(), "src", "data", "alignments")
        files = os.listdir(dir_path)
        alignments = {}
        for file in files:
            file_path = os.path.join(os.getcwd(), "src", "data", "alignments", file)
            gene_name = file.split("_")[self.GENE_NAME_IDX]
            with open(file_path, "r") as f:
                if not alignments.get(gene_name):
                    alignments[gene_name] = {}
                org = ""  # keep track of the last seen organism
                for line in f.readlines():
                    if line.startswith(">"):
                        org = line.replace(">", "").strip()
                        if not alignments.get(gene_name).get(org):
                            alignments[gene_name][org] = ""
                    else:  # we're on an alignment line
                        alignments[gene_name][org] = alignments[gene_name][org] + line.strip()
        return alignments

    def _get_significant_orgs(self):
        """ Reads in the significant organisms selected for this study """
        orgs = {}
        path = os.path.join(os.getcwd(), "src", "data", "phylogeny", "significant_organisms.txt")
        with open(path, "r") as f:
            f.readline()  # don't care about the top line
            for line in f.readlines():
                org_name = line.split(",")[self.ORG_NAME_IDX]
                org = "_".join(org_name.lower().split())
                orgs[org] = 1
        return orgs

    def _setup_R_access(self):
        """ Sets up the necessary R components such as libs and functions """
        utils = packages.importr('utils')
        # select a mirror for R packages
        utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
        # R package names
        packnames = ['ape']
        names_to_install = [x for x in packnames if not packages.isinstalled(x)]
        if len(names_to_install) > 0:
            utils.install_packages(vectors.StrVector(names_to_install))
        robjects.r(
            '''
            # function to compute the DN/DS ratio of two sequences
            library(ape)
            get_dnds = function(seq1, seq2) {
                mtx = c(strsplit(seq1, ""), strsplit(seq2, ""))
                bin = as.DNAbin(mtx)
                dnds(bin)
            }
            '''
        )
        return robjects.r['get_dnds']

    def _get_dnds_scores(self):
        """ Computes and returns a dictionary of all the dnds scores for each pair of organisms for a gene """
        final = {}
        # dnds < 1 = negative selection, dnds == 1 = neutral, dnds > 1 = positive selection
        for org1 in self.significant_orgs.keys():
            for org2 in self.significant_orgs.keys():
                for gene in self.alignments.keys():
                    if not final.get(org1):
                        final[org1] = {}
                    if not final.get(org1).get(org2):
                        final[org1][org2] = {}
                    if org1 == org2:
                        final[org1][org2][gene] = 1  # neutral for the same organism
                    else:
                        alignment = self.alignments.get(gene)
                        seq1 = alignment.get(org1)
                        seq2 = alignment.get(org2)
                        if not seq1 or not seq2:
                            final[org1][org2][gene] = 1
                            # in this case, one of the organisms does not have a gene
                            # assume the score is 1 for neutrality
                            continue
                        seq1_clean, seq2_clean = self._clean_sequence(seq1, seq2, check_chars=True)
                        # apparently, it can happen that we get non-unique sequences after cleaning
                        if seq1_clean == seq2_clean:
                            final[org1][org2][gene] = 1
                        else:
                            try:
                                score = self.get_dnds(seq1_clean, seq2_clean)[0]  # results come back as vectors
                                # placing a 1 for neutrality when we get NaNs, based on docs, the sequences should be
                                # highly divergent but don't know for sure, safer to place in neutral, might change
                                # mind though...
                                final[org1][org2][gene] = round(score, 2) if not np.isnan(score) else 1
                            except Exception as e:
                                print("CAUGHT EXCEPTION -- {}".format(e))
                                # exceptions can be thrown if two sequences are very divergent, assume a 1 for neutral
                                final[org1][org2][gene] = 1
        return final

    def _clean_sequence(self, seq1, seq2, check_chars=False):
        """
        Cleans the given pair of sequences such that dnds can be computed
        :param str seq1: first sequence to clean
        :param str seq2: second sequence to clean
        """
        assert (len(seq1) == len(seq2))
        l_seq1 = list(seq1)
        l_seq2 = list(seq2)
        for i in range(len(l_seq1)):  # or whatever
            gap = l_seq1[i] == "-" or l_seq2[i] == "-"
            if gap:
                l_seq1[i], l_seq2[i] = "@", "@"
            elif check_chars:
                # there are sequences that contain non-IUPAC char, which cannot be used for dnds evaluation
                illegal = self._is_illegal_char(l_seq1[i]) or self._is_illegal_char(l_seq2[i])
                if illegal:
                    l_seq1[i], l_seq2[i] = "@", "@"
        seq1, seq2 = "".join(l_seq1).replace("@", ""), "".join(l_seq2).replace("@", "")

        trim_len = min(len(seq1), len(seq2))
        seq1_trimmed = seq1[:trim_len]
        seq2_trimmed = seq2[:trim_len]

        seq1_remainder = len(seq1_trimmed) % 3
        if not seq1_remainder == 0:
            seq1_trimmed = seq1_trimmed[:len(seq1_trimmed) - seq1_remainder]
        seq2_remainder = len(seq2_trimmed) % 3
        if not seq2_remainder == 0:
            seq2_trimmed = seq2_trimmed[:len(seq2_trimmed) - seq2_remainder]
        assert (len(seq1_trimmed) == len(seq2_trimmed))
        return seq1_trimmed, seq2_trimmed

    def _is_illegal_char(self, c):
        """ Check if a character is illegal - non-IUPAC nucleotides """
        if c in self.UNACCEPTED_CHARS:
            return True
        return False
