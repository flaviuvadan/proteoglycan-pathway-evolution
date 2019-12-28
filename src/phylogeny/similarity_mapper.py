import os


class Mapper:
    """ Responsible for finding similarities between phylogenetic trees and creating visualizations for the results """

    GENE_IDX = 0
    ORGANISMS_IDX = 1

    SIMILARITY_THRESHOLD = 0.85

    def __init__(self):
        self.genes_organisms = self._load_tree_data()

    def _load_tree_data(self):
        """
        Builds and returns a dictionary mapping genes to the organisms they appear in
        :return: dictionary of gene organisms pairs
        """
        genes_organisms = {}
        with open(os.path.join(os.getcwd(), "src", "data", "genes", "genes_organisms.txt"), "r") as f:
            f.readline()  # don't care about the first line
            for line in f.readlines():
                split_line = line.split(",")
                gene, organisms = split_line[self.GENE_IDX], split_line[self.ORGANISMS_IDX]
                split_organisms = organisms.split("/")
                genes_organisms[gene] = {}
                for org in split_organisms:
                    genes_organisms[gene][org] = 1
        return genes_organisms

    def _find_similar_tree(self, gene_org, gene_orgs):
        """
        Finds the highest similarity tree for the given tree
        :param dict gene_org: a a gene
        :param list gene_orgs: list of pairs of (gene, orgs)
        :return: a tuple of (similar tree key, similar tree orgs, similarity score)
        """
        gene_org_key, gene_org_orgs = gene_org[self.GENE_IDX], gene_org[self.ORGANISMS_IDX]
        max_similarity = 0
        similar_tree = None
        for go in gene_orgs:
            go_orgs = go[self.ORGANISMS_IDX]
            intersection = set(gene_org_orgs).intersection(go_orgs)
            union = set(gene_org_orgs).union(set(go_orgs))
            similarity = len(intersection) / len(union)  # Jaccard's index
            if similarity > max_similarity:
                max_similarity = similarity
                similar_tree = (go[self.GENE_IDX], go[self.ORGANISMS_IDX])
        # this will cause the queue to remove two trees since the score is 0, which allows it to decrease in size
        if not similar_tree:
            return None, None, 0
        return similar_tree[self.GENE_IDX], similar_tree[self.ORGANISMS_IDX], max_similarity

    def _create_similarity_trees(self):
        """
        Constructs and returns a list of combined gene organisms based on how similar the orgs are
        :return: a list of new (gene, organism)s
        """
        new_trees = []
        queue = []
        for k, v in self.genes_organisms.items():
            queue.append((k, v))
        while len(queue) > 0:
            pair = queue.pop(0)  # element at the head of the queue
            gene, orgs = pair[self.GENE_IDX], pair[self.ORGANISMS_IDX]
            similar_key, similar, score = self._find_similar_tree(pair, queue)
            if score >= self.SIMILARITY_THRESHOLD:
                new_key = os.path.join(gene, similar_key)
                new_orgs = set(orgs).intersection(set(similar))
                try:
                    to_remove = (similar_key, similar)
                    queue.remove(to_remove)
                except ValueError:
                    # this is weird if it happens
                    print("something is funny")
                    continue
                queue.append((new_key, new_orgs))
            else:
                new_trees.append((gene, set(orgs)))
        return new_trees

    def visualize_trees(self):
        similarity_trees = self._create_similarity_trees()
        for st in similarity_trees:
            print("{}\n{}".format(st[0], st[1]))


if __name__ == "__main__":
    mapper = Mapper()
    mapper.visualize_trees()
