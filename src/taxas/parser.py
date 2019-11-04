import csv
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn


class Parser:
    """ The parser is responsible for parsing the organism gene frequencies and the gene collections for plotting """

    GENE_NAME_IDX = 0
    GENE_ID_IDX = 1
    GENE_IDX = 2

    def __init__(self):
        """ Constructor """
        self.genes, self.gene_indices = self.load_genes()
        self.gene_freqs = self._parse_gene_freqs()
        self.org_genes = self._parse_organisms_genes()

    def load_genes(self):
        """ Reads in the genes """
        genes = []
        genes_indices = {}
        with open(os.path.join(os.getcwd(), "src", "data", "genes", "genes.txt"), 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for i, gene in enumerate(csv_reader):
                # record the i'th index because we want to keep track of the location of the gene in the binary vector
                # see _parse_organisms_genes()
                genes_indices[gene[self.GENE_NAME_IDX]] = i
                genes.append((gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX]))
        return genes, genes_indices

    def _parse_gene_freqs(self):
        """ Parses the gene frequencies file """
        frequencies = {}
        with open(os.path.join(os.getcwd(), "src", "data", "genes", "gene_frequencies.txt"), "r") as freqs:
            freqs.readline()  # omit the first line
            for line in freqs.readlines():
                split_line = line.split(",")
                org_name, gene_freq = split_line[0], int(split_line[1].replace("\n", ""))
                frequencies[org_name] = gene_freq
        return frequencies

    def create_gene_freq_plot(self):
        """ Creates a plot of the gene frequencies """
        plt.figure(figsize=(20, 20))
        plt.scatter(self.gene_freqs.values(), self.gene_freqs.keys(), marker='.')
        plt.title("Organisms gene frequencies", fontsize=40)
        plt.axis([21, max(self.gene_freqs.values()) + 2, 0, len(self.gene_freqs.keys())])
        plt.yticks(range(len(list(self.gene_freqs.keys()))), list(self.gene_freqs.keys()))
        plt.tick_params(axis='y', labelsize=6)
        plt.savefig("gene_frequencies.pdf", format="pdf", quality=95)
        plt.close()

    def _get_genes_vector(self):
        """ Builds and returns a gene binary vector indicating the presence or absence of a gene. The indices
        correspond to the way genes are loaded from genes.txt """
        return [0] * len(self.genes)

    def _parse_organisms_genes(self):
        """ Parses the organisms genes file """
        org_genes = {}
        # add each org to dict, each pointing to another dict that's a binary gene representation
        with open(os.path.join(os.getcwd(), "src", "data", "genes", "organisms_genes.txt"), "r") as orgs:
            orgs.readline()  # omit the first line
            for line in orgs.readlines():
                split_line = line.split(",")
                organism, genes = split_line[0], split_line[1].strip("\n").split(" ")
                org_genes[organism] = self._get_genes_vector()
                for gene in genes:
                    if not gene:
                        continue
                    try:
                        org_genes[organism][self.gene_indices.get(gene)] = 1
                    except Exception as e:
                        print("Caught exception: {}, gene: {}".format(e, gene))
                        break
        return org_genes

    def _create_organisms_genes_csv(self):
        """ Creates a csv of the generated organisms gene vectors """
        with open("organisms_genes_vectors.txt", "w") as f:
            title = "Organism,"
            total_gene_indices = len(self.gene_indices.keys())
            for i, gene in enumerate(self.gene_indices.keys()):
                title = title + ("{},".format(gene) if i != total_gene_indices - 1 else "{}".format(gene))
            f.write(title + "\n")
            for org, vector in self.org_genes.items():
                line = "{},".format(org)
                vector_len = len(vector)
                for i, v in enumerate(vector):
                    line = line + "{},".format(v) if i != vector_len - 1 else line + "{}".format(v)
                line = line + "\n"
                f.write(line)

    def create_organisms_genes_matrix(self):
        """ Creates a plot that is a binary matrix of all the genes vs. organisms """
        # easier to go about this with this approach than fight with matplotlib
        path = os.path.join(os.getcwd(), "src", "data", "genes", "organisms_genes_vectors.txt")
        df = pd.read_csv(path)
        y_values = df.Organism
        x_values = df.columns[1:]
        df = df.drop(columns=["Organism"])

        fig_text = "A matrix representation of the presence of each gene in every studied organism"

        plt.figure(figsize=(50, 50))
        plt.title("Organisms gene presence", fontsize=50)
        plt.figtext(0.5, 0.01, fig_text, wrap=True, horizontalalignment="center", fontsize=40)
        plt.xticks(range(len(x_values)), x_values)
        plt.yticks(range(len(y_values)), y_values)

        sn.heatmap(df, cbar=False, xticklabels=x_values, yticklabels=y_values, cmap="Blues")
        plt.savefig("gene_presence_vectors.pdf")


if __name__ == "__main__":
    parser = Parser()
    parser.create_gene_freq_plot()
    parser.create_organisms_genes_matrix()
