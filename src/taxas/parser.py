import csv
import os


class Parser:
    """ The parser is responsible for parsing the organism gene frequencies and the gene collections for plotting """

    GENE_ID_IDX = 1
    GENE_NAME_IDX = 0

    def __init__(self, genes_file_path):
        """ Constructor """
        self.genes = self.load_genes(genes_file_path)

    def load_genes(self, genes_file_path):
        """ Reads in the genes """
        with open(genes_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield (gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])

    def parse(self):
        """ Parses the two files - class doc """
        pass

    def _parse_gene_freqs(self):
        """ Parses the gene frequencies file """
        pass

    def _parse_organisms_genes(self):
        """ Parses the organisms genes file """
        pass

    def


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    parser = Parser(gene_file_path)
