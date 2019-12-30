import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import packages, vectors
import seaborn as sn


class Visualizer:
    """ Responsible for creating dn/ds visualizations for all the studied genes """

    GENE_NAME_IDX = 0

    ORG_CLASS_IDX = 0
    ORG_NAME_IDX = 1

    UNACCEPTED_CHARS = ["W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]

    def __init__(self):
        """ Constructor """
        self.alignments = self._parse_organisms_alignments()
        self.significant_orgs = self._get_significant_orgs()
        self.get_dnds = self._setup_R_access()
        self.scores = self._get_dnds_scores()
        self._create_heatmap_csv()

    def _create_heatmap_csv(self):
        """ Creates a heatmap CSV file for visualizing the results """
        for sig_org in self.scores.keys():
            sig_org_file_path = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(sig_org))
            with open(sig_org_file_path, "w") as sof:
                header = "organism,"
                # each org will have itself in the dnds dictionary
                # we only want this to build the header once
                for k in self.scores[sig_org][sig_org].keys():
                    header += "{},".format(k)
                header = header[:len(header) - 1]  # remove the last ","
                header += "\n"
                sof.write(header)
                for org in self.scores[sig_org].keys():
                    if sig_org == org:  # no point in adding the same one
                        continue
                    vals = ",".join([str(x) for x in self.scores[sig_org][org].values()])
                    line = "{},{}\n".format(org, vals)
                    sof.write(line)

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
                        seq1_clean, seq2_clean = self._clean_sequence(seq1, seq2,
                                                                      check_chars=True)
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
                            except Exception:
                                final[org1][org2][gene] = 0
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
            illegal = l_seq1[i] == "-" or l_seq2[i] == "-"
            if check_chars:
                # there are sequences that contain non-IUPAC char, which cannot be used for dnds evaluation
                illegal = self._is_illegal_char(l_seq1[i]) or self._is_illegal_char(l_seq2[i])
                if illegal:
                    l_seq1[i], l_seq2[i] = "@", "@"
            elif illegal:
                # mark characters for removal, this is necessary as strings are immutable
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

    def visualize(self):
        """ Creates the 2D dn/ds plots for each gene """
        for sig_org in self.scores.keys():
            sig_org_csv = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(sig_org))
            df = pd.read_csv(sig_org_csv)
            x_vals = df.columns[1:]
            y_vals = df.organism
            y_vals_clean = [" ".join(y.split("_")).capitalize() for y in y_vals]
            df = df.drop(columns="organism")
            plt.figure(figsize=(20, 20))
            plt.title("{} DN/DS ratios".format(" ".join(sig_org.split("_")).capitalize()), fontsize=30)
            plt.xticks(range(len(x_vals)), x_vals,
                       rotation=45)
            plt.yticks(range(len(y_vals_clean)), y_vals_clean,
                       rotation=45)
            sn.heatmap(df,
                       cbar=True,
                       xticklabels=x_vals,
                       yticklabels=y_vals_clean,
                       cmap="YlGnBu",
                       fmt="f",
                       linewidths=0.5)
            save_fig_path = os.path.join(os.getcwd(), "src", "data", "visualizations", "dnds", "{}.pdf".format(sig_org))
            plt.savefig(save_fig_path,
                        quality=95,
                        bbox_inches='tight')
            plt.clf()

    def _contains_illegal_chars(self, seq):
        """ Check whether a given sequence contains non-ATCG letters (IUPAC accepted but cannot compute dn/ds with
        those) """
        for u in self.UNACCEPTED_CHARS:
            if u in seq:
                return True
        return False

    def _is_illegal_char(self, c):
        """ Check if a character is illegal - non-IUPAC nucleotides """
        if c in self.UNACCEPTED_CHARS:
            return True
        return False


if __name__ == "__main__":
    visualizer = Visualizer()
    visualizer.visualize()
