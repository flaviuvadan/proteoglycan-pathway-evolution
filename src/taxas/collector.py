import csv
import os


class Collector:
    """ Collector is responsible for collecting taxonomic information for all the unique organisms in data/organisms"""

    GENE_ID_IDX = 1

    def __init__(self, gene_file_path):
        """
        Constructor
        :param gene_file_path: name of the file containing gene names
        """
        self.gene_file_path = gene_file_path
        self.gene_ids = self.load_genes()
        self.organisms = self.load_organisms()

    def load_genes(self):
        """ Loads the genes into the class gene_ids """
        with open(self.gene_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield gene[self.GENE_ID_IDX]

    def get_organisms_file_path(self, gene_id):
        """ Builds the file path to the organisms file of the given gene ID """
        return os.path.join(os.pardir, "data", "organisms", "{}.txt".format(gene_id))

    def load_organisms(self):
        """ Responsible for loading all the organisms curated for the list of genes of interest """
        organisms = {}
        for gene_id in self.gene_ids:
            org_file_path = self.get_organisms_file_path(gene_id)
            with open(org_file_path, "r") as orgs:
                org = orgs.read().splitlines()
                # we only care about unique organisms
                for o in org:
                    if not organisms.get(o):
                        organisms[o] = 1
                    else:
                        organisms[o] = organisms.get(o) + 1
        return organisms

    def collect(self):
        """ Collects taxa information from GBIF (https://www.gbif.org/developer/species) """
        print("would collect")
        pass


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    collector = Collector(gene_file_path)
    collector.collect()
