import os

import ete3 as ete


class Mapper:
    """ Responsible for finding similarities between phylogenetic trees and creating visualizations for the results """

    GENE_IDX = 0
    ORGANISMS_IDX = 1

    GENUS_IDX = 0
    SPECIES_IDX = 1

    SIMILARITY_THRESHOLD = 0.80

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

    def _visualise_similarity_tree(self, newick_info, gene, organisms):
        """
        Creates the radial gene tree of each gene
        :param newick_info - newick-formatted tree representation
        :param gene - gene name
        :param organisms - list of organisms that contain the given gene
        """
        t = ete.Tree(newick_info, format=1)

        ts = ete.TreeStyle()
        ts.show_leaf_name = True
        ts.mode = 'c'
        ts.arc_span = 360
        ts.force_topology = True
        ts.title.add_face(ete.TextFace(gene, fsize=150, bold=True), column=0)
        ts.show_leaf_name = False
        ts.show_scale = False

        ns = ete.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1

        for node in t.traverse():
            node.set_style(ns)
            if node.is_leaf():
                split_node_name = node.name.strip("'").split()
                genus_species = " ".join([split_node_name[self.GENUS_IDX], split_node_name[self.SPECIES_IDX]])
                found_in_orgs = False
                for org in organisms:
                    if genus_species.lower() in org.lower():
                        found_in_orgs = True
                        break
                if found_in_orgs:
                    node_face = ete.TextFace(node.name.strip("'"), fsize=60, penwidth=10)
                    node.add_face(node_face, 1)
                else:
                    node_face = ete.TextFace(node.name.strip("'"), fsize=55, penwidth=10, fgcolor="blue")
                    node.add_face(node_face, 1)
        file_name = "_".join(gene.split("/"))
        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "similarity",
                                   "{}.pdf".format(file_name))
        t.render(destination, tree_style=ts)

    def visualize_trees(self):
        """ Creates radial trees for visualising similar, combined, trees """
        similarity_trees = self._create_similarity_trees()
        path_to_newick = os.path.join(os.getcwd(), "src", "data", "phylogeny", "all_orgs_phyliptree.phy")
        with open(path_to_newick, "r") as f:
            newick_info = f.read()
            for idx, tree in enumerate(similarity_trees):
                gene, orgs = tree[self.GENE_IDX], tree[self.ORGANISMS_IDX]
                self._visualise_similarity_tree(newick_info, gene, orgs)


if __name__ == "__main__":
    mapper = Mapper()
    mapper.visualize_trees()
