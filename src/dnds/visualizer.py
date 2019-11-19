import os

import dnds
import pandas as pd


class Visualizer:
    """ Responsible for creating dn/ds visualizations for all the studied genes """

    GENE_NAME_IDX = 0

    def __init__(self):
        """ Constructor """
        self.alignments = self._parse_organisms_alignments()
        self.dnds_dataframes, self.error_log = self._get_gene_dnds_dataframes()

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

    def _get_gene_dnds_dataframes(self):
        """ Builds and returns a list of dataframes that contain the dn/ds scores of each gene """
        dfs = {}
        went_wrong = []
        for gene_key in self.alignments.keys():
            gene = self.alignments.get(gene_key)
            dfs[gene_key] = pd.DataFrame(index=list(gene.keys()), columns=list(gene.keys()))
            for org_1 in gene.keys():
                for org_2 in gene.keys():
                    org_1_seq, org_2_seq = gene.get(org_1), gene.get(org_2)
                    if org_1_seq == org_2_seq:
                        orgs_dnds = 0
                    elif self._contains_illegal_letters(org_1_seq) or self._contains_illegal_letters(org_2_seq):
                        msg = "gene: {}, org1: {}, org2: {}, exception: {}".format(gene_key, org_1, org_2,
                                                                                   "illegal characters")
                        orgs_dnds = 0
                    else:
                        org_1_seq_list = list(org_1_seq)
                        org_2_seq_list = list(org_2_seq)
                        # expect these to be equal because of the MSA
                        assert (len(org_1_seq_list) == len(org_2_seq_list))
                        for i in range(len(org_1_seq_list)):
                            if org_1_seq_list[i] == "-" or org_2_seq_list[i] == "-":
                                org_1_seq_list[i], org_2_seq_list[i] = "*", "*"
                        # these should have the same modulo remainder
                        org_1_seq = "".join(org_1_seq_list).replace("*", "")
                        org_2_seq = "".join(org_2_seq_list).replace("*", "")
                        if len(org_1_seq) % 3 != 0:
                            org_1_seq = org_1_seq[:-(len(org_1_seq) % 3)]
                        if len(org_2_seq) % 3 != 0:
                            org_2_seq = org_2_seq[:-(len(org_2_seq) % 3)]
                        assert (len(org_1_seq) == len(org_2_seq))
                        try:
                            orgs_dnds = dnds.dnds(org_1_seq, org_2_seq)
                        # a bunch of shit can go wrong..
                        except Exception as e:
                            msg = "gene: {}, org1: {}, org2: {}, exception: {}".format(gene_key, org_1, org_2, e)
                            went_wrong.append(msg)
                            orgs_dnds = 0
                    dfs[gene_key][org_1][org_2] = orgs_dnds
            break
        return dfs, went_wrong

    def _contains_illegal_letters(self, seq):
        """ Check whether a given sequence contains non-ATCG letters (IUPAC accepted but cannot compute dn/ds with
        those) """
        unaccepted = ["W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]
        for u in unaccepted:
            if u in seq:
                return True
        return False

    def visualize(self):
        """ Creates the 2D dn/ds plots for each gene """
        for k1 in self.alignments.keys():
            for k2 in self.alignments.get(k1).keys():
                print("{} -- {}".format(k1, k2))
                print(self.alignments.get(k1).get(k2))
            break


if __name__ == "__main__":
    visualizer = Visualizer()
    for err in visualizer.error_log:
        print(err)
    # visualizer.visualize()
