import os

import dnds
import pandas as pd


class Visualizer:
    """ Responsible for creating dn/ds visualizations for all the studied genes """

    GENE_NAME_IDX = 0

    def __init__(self):
        """ Constructor """
        self.alignments = self._parse_organisms_alignments()
        self.dnds_dataframes = self._get_gene_dnds_dataframes()

    def _parse_organisms_alignments(self):
        """ Reads in the organisms alignments """
        dir_path = os.path.join(os.getcwd(), "src", "data", "organisms")
        files = os.listdir(dir_path)
        alignments = {}
        for file in files:
            file_path = os.path.join(os.getcwd(), "src", "data", "organisms", file)
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
                        alignments[gene_name][org] = line.strip()
        return alignments

    def _get_gene_dnds_dataframes(self):
        """ Builds and returns a list of dataframes that contain the dn/ds scores of each gene """
        dfs = []
        for gene_key in self.alignments:
            for org_1 in self.alignments.get(gene_key).keys():
                for org_2 in self.alignments.get(gene_key).keys():
                    pass
        return dfs

    def _compute_dnds(self, seq1, seq2):
        """ Compute the dnds score between two sequences """
        return dnds.dnds(seq1, seq2)

    def visualize(self):
        """ Creates the 2D dn/ds plots for each gene """
        for k1 in self.alignments.keys():
            for k2 in self.alignments.get(k1).keys():
                print("{} -- {}".format(k1, k2))
                print(self.alignments.get(k1).get(k2))
            break


if __name__ == "__main__":
    visualizer = Visualizer()
    visualizer.visualize()
