import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn


class Mapper:
    """ Mapper is responsible for holding logic related to building maps/heatmaps of genetic data and
    parsing/organizing data for clustering """

    def __init__(self):
        """ Constructor """
        self.base_path = os.getcwd()
        self.gene_vectors = self._load_gene_vectors()

    def _load_gene_vectors(self):
        """ Loads the organisms gene vectors from a data file """
        return pd.read_csv(self._get_gene_vectors_file_path())

    def _get_gene_vectors_file_path(self):
        """ Builds and returns the file path to the organisms gene vectors """
        return os.path.join(self.base_path, "data", "organisms_genes_vectors.txt")

    def build_hierchical_heatmap(self):
        """ Builds a hierarchical heatmap of the gene vectors """
        y_values = self.gene_vectors.Organism
        x_values = self.gene_vectors.columns[1:]
        df = self.gene_vectors.drop(columns="Organism")
        fig_text = "A clustermap of organisms vs. genes representing the Euclidean between the gene vectors"
        # plt.figure(figsize=(50, 50))
        plt.title("Organisms gene presence", fontsize=50)
        plt.figtext(0.5, 0.01, fig_text, wrap=True, horizontalalignment="center", fontsize=40)
        plt.xticks(range(len(x_values)), x_values)
        plt.yticks(range(len(y_values)), y_values)
        fig = sn.clustermap(df, cbar=False, xticklabels=x_values, yticklabels=y_values, cmap=["White", "Blue"], figsize=(50,50))
        fig.cax.set_visible(False)
        fig.savefig("genes_clustermap.pdf")


if __name__ == "__main__":
    mapper = Mapper()
    mapper.build_hierchical_heatmap()
