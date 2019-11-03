import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn


class Visualizer:
    """ Responsible for creating a simple visualization of the gene functions """

    def __init__(self):
        """ Constructor """
        self.genes = self._read_genes()

    def _read_genes(self):
        """ Reads in the genes CSV file to create a dataframe """
        return pd.read_csv("gene_functions.txt")

    def visualize(self):
        """ Creates the grid representing gene functions """
        x_vals = self.genes.columns[1:]
        y_vals = self.genes.Gene
        df = self.genes.drop(columns="Gene")

        plt.figure(figsize=(20, 20))
        plt.title("Matrix representation of the gene functions", fontsize=30)
        plt.xticks(range(len(x_vals)), x_vals, rotation="45")
        plt.yticks(range(len(y_vals)), y_vals)
        sn.heatmap(df, cbar=False, xticklabels=x_vals, yticklabels=y_vals, cmap=["White", "Blue"])
        plt.savefig("gene_functions.pdf", quality=95)


if __name__ == "__main__":
    viz = Visualizer()
    viz.visualize()
