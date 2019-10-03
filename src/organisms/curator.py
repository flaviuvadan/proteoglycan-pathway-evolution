import csv
import json
import os

from src import exceptions


class Curator:
    """ Curator is responsible for loading the 51 genes of interest, reading their orthologs, and parsing their
    Ensembl results with the intention to create file with unique organisms for those genes """

    GENE_ID_IDX = 1

    def __init__(self, gene_file_path):
        """
        Constructor
        :param gene_file_path: name of the file containing gene names
        """
        self.gene_file_path = gene_file_path
        self.gene_ids = self.load_genes()

    def load_genes(self):
        """ Loads the genes into the class gene_ids """
        with open(self.gene_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield gene[self.GENE_ID_IDX]

    def build_gene_file_path(self, gene_id):
        """
        Builds the path to a gene orthologs file
        :param gene_id: Ensembl ID of the gene
        :return: path to file
        """
        return os.path.join(os.pardir, "data", "orthologs", "{}.txt".format(gene_id))

    def curate(self):
        """ Curates organisms """
        for gene_id in self.gene_ids:
            organisms = {}
            with open(self.build_gene_file_path(gene_id), "r") as orthologs:
                loaded = json.load(orthologs)
                data = loaded.get('data')
                if not data:
                    raise exceptions.EmptyOrthologData("gene {} orthologs not found".format(gene_id))
                homologies = data[0].get('homologies')
                if not homologies:
                    raise exceptions.EmptyOrthologData("gene {} ortholog homologies not found".format(gene_id))
                for homology in homologies:
                    source = homology.get('source')
                    if not source:
                        raise exceptions.EmptyHomologyInformation("gene {} ortholog has no source".format(gene_id))
                    target = homology.get('target')
                    if not target:
                        raise exceptions.EmptyHomologyInformation("gene {} ortholog has no target".format(gene_id))
                    source_species = source.get('species')
                    if not source_species:
                        raise exceptions.EmptyHomologyInformation(
                            "gene {} ortholog has not source target".format(gene_id))
                    target_species = target.get('species')
                    if not target_species:
                        raise exceptions.EmptyHomologyInformation(
                            "gene {} ortholog has not source species".format(gene_id))
                    # we only need to get all the unique organisms
                    organisms[source_species] = 1
                    organisms[target_species] = 1
            self.create_organisms_file(gene_id, organisms)
            break

    def get_organisms_file_path(self, gene_id):
        """ Builds the file path to the organisms file of the given gene ID """
        return os.path.join(os.pardir, "data", "organisms", "{}.txt".format(gene_id))

    def create_organisms_file(self, gene_id, organisms):
        """
        Creates the organisms files associated with the orthologs of a gene
        :param str gene_id: Ensembl ID of the gene
        :param dict organisms: dictionary of organisms keyed on species
        """
        organisms_file_path = self.get_organisms_file_path(gene_id)
        with open(organisms_file_path, 'w') as out:
            for species in organisms.keys():
                out.write(species.replace("_", " ") + "\n")


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    curator = Curator(gene_file_path)
    curator.curate()
