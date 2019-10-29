import csv
import os

import matplotlib.pyplot as plt


class Parser:
    """ The parser is responsible for parsing the organism gene frequencies and the gene collections for plotting """

    GENE_ID_IDX = 1
    GENE_NAME_IDX = 0

    def __init__(self, genes_file_path):
        """ Constructor """
        self.genes = self.load_genes(genes_file_path)
        self.gene_freqs = self._parse_gene_freqs()

    def load_genes(self, genes_file_path):
        """ Reads in the genes """
        with open(genes_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield (gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])

    def _parse_gene_freqs(self):
        """ Parses the gene frequencies file """
        frequencies = {}
        with open("gene_frequencies.txt", "r") as freqs:
            freqs.readline()  # omit the first line
            for line in freqs.readlines():
                split_line = line.split(",")
                org_name, gene_freq = split_line[0], int(split_line[1].replace("\n", ""))
                frequencies[org_name] = gene_freq
        return frequencies

    def create_gene_freq_plot(self):
        """ Creates a plot of the gene frequencies """
        plt.figure(figsize=(18, 18))
        plt.scatter(self.gene_freqs.keys(), self.gene_freqs.values(), marker='.')
        plt.title("Organisms gene frequencies")
        plt.axis([0, len(self.gene_freqs.keys()), 20, max(self.gene_freqs.values()) + 3])
        plt.xticks(list(self.gene_freqs.keys()), rotation='vertical')
        plt.tick_params(axis='x', which='major', labelsize=5)
        plt.savefig("gene_freqs.pdf", format="pdf", quality=95)
        plt.close()

    def _parse_organisms_genes(self):
        """ Parses the organisms genes file """
        pass


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    parser = Parser(gene_file_path)
    parser.create_gene_freq_plot()
