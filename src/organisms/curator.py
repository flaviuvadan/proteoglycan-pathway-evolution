import csv
import json
import os

from src import exceptions


class Curator:
    """ Curator is responsible for loading the 51 genes of interest, reading their orthologs, and parsing their
    Ensembl results with the intention to create file with unique organisms for those genes """

    GENE_NAME_IDX = 0
    GENE_ID_IDX = 1

    def __init__(self, gene_file_path):
        """
        Constructor
        :param gene_file_path: name of the file containing gene names
        """
        self.gene_file_path = gene_file_path
        self.genes = self.load_genes()

    def load_genes(self):
        """ Loads the genes into the class gene_ids """
        with open(self.gene_file_path, 'r') as gene_file:
            csv_reader = csv.reader(gene_file, delimiter=',')
            for gene in csv_reader:
                yield (gene[self.GENE_NAME_IDX], gene[self.GENE_ID_IDX])

    def build_gene_file_path(self, gene_id):
        """
        Builds the path to a gene orthologs file
        :param gene_id: Ensembl ID of the gene
        :return: path to file
        """
        return os.path.join(os.pardir, "data", "orthologs", "{}.txt".format(gene_id))

    def curate(self):
        """ Curates organisms """
        for gene in self.genes:
            gene_id = gene[self.GENE_ID_IDX]
            gene_name = gene[self.GENE_NAME_IDX]
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
                    source_species, source_seq, target_species, target_seq = self.parse_homologies(gene_id, homology)
                    # want to keep seqs around for multiple sequence alignment
                    organisms[source_species] = source_seq
                    organisms[target_species] = target_seq
            self.create_organisms_file(gene_name, gene_id, organisms)

    def parse_homologies(self, gene_id, homology):
        """
        Parses the homologies dictionary of the orthologs of a gene to get the source and target species, along with
        their sequences
        :param str gene_id: Ensembl ID of the gene
        :param dict homology: information of a gene
        :return: source species, target species
        """
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
        source_seq = source.get('align_seq')
        if not source_seq:
            raise exceptions.EmptySquence("gene {} source seq not found".format(gene_id))
        target_seq = target.get('align_seq')
        if not target_seq:
            raise exceptions.EmptySquence("gene {} target seq not found".format(gene_id))
        return source_species, source_seq, target_species, target_seq

    def get_organisms_file_path(self, gene_name, gene_id):
        """ Builds the file path to the organisms file of the given gene name and gene ID """
        return os.path.join(os.pardir, "data", "organisms", "{}_{}.txt".format(gene_name, gene_id))

    def create_organisms_file(self, gene_name, gene_id, organisms):
        """
        Creates the organisms files associated with the orthologs of a gene
        :param str gene_name: name of the gene being processed
        :param str gene_id: Ensembl ID of the gene
        :param dict organisms: dictionary of organisms keyed on species
        """
        organisms_file_path = self.get_organisms_file_path(gene_name, gene_id)
        with open(organisms_file_path, 'w') as out:
            for species, sequence in organisms.items():
                out.write("{}:{}\n".format(species, sequence.replace("-", "")))


if __name__ == "__main__":
    gene_file_path = os.path.join(os.pardir, "data", "genes.txt")
    curator = Curator(gene_file_path)
    curator.curate()
