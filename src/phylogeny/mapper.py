import os

import ete3
from ete3 import Tree, TreeStyle


class Mapper:
    """ Mapper is responsible for creating phylogenetic trees for all the genes based on a global phylogenetic tree
    created based on taxonomic information of all the organisms used in this study """

    GENE_IDX = 0
    ORGANISMS_IDX = 1

    def __init__(self):
        """ Constructor """
        self.genes_organisms = self._read_genes_organisms()

    def _read_genes_organisms(self):
        """ Reads in the gene organisms file and creates a lookup dictionary """
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

    def create_global_trees(self):
        """ Creates the main global tree """
        with open(os.path.join(os.getcwd(), "src", "data", "phylogeny", "phyliptree.phy"), "r") as f:
            newick_info = f.read()
            self._create_radial_tree(newick_info)
            self._create_standard_tree(newick_info)

    def _create_radial_tree(self, newick_info):
        """
        Creates the radial version of the global tree
        :param newick_info - newick-formatted tree representation
        """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.mode = 'c'
        ts.arc_span = 360
        ts.force_topology = True
        ts.legend = None

        ns = ete3.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1
        for node in t.traverse():
            node.set_style(ns)
            if node.is_leaf():
                node_face = ete3.TextFace(node.name.strip("'"), fsize=60, penwidth=10)
                node.add_face(node_face, 1)

        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "radial", "global.pdf")
        t.render(destination, tree_style=ts)

    def _create_standard_tree(self, newick_info):
        """
        Creates the standard version of the global tree
        :param newick_info - newick-formatted tree representation
        """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.force_topology = True
        ts.legend = None
        ts.rotation = 90

        ns = ete3.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1
        for node in t.traverse():
            node.set_style(ns)
            if node.is_leaf():
                node_face = ete3.TextFace(node.name.strip("'"), fsize=60, penwidth=10)
                node.add_face(node_face, 1)

        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "standard", "global.pdf")
        t.render(destination, tree_style=ts)

    def create_gene_trees(self):
        """ Creates the 51 trees, one per gene """
        with open(os.path.join(os.getcwd(), "src", "data", "phylogeny", "phyliptree.phy"), "r") as f:
            newick_info = f.read()
            for gene in self.genes_organisms.keys():
                self._create_radial_gene_tree(newick_info, gene, self.genes_organisms.get(gene).keys())
                self._create_standard_gene_tree(newick_info, gene, self.genes_organisms.get(gene).keys())

    def _create_radial_gene_tree(self, newick_info, gene, organisms):
        """
        Creates the radial gene tree of each gene
        :param newick_info - newick-formatted tree representation
        :param gene - gene name
        :param organisms - list of organisms that contain the given gene
        """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = 'c'
        ts.arc_span = 360
        ts.force_topology = True
        ts.title.add_face(ete3.TextFace(gene, fsize=120, bold=True), column=0)
        ts.show_leaf_name = False
        ts.show_scale = False

        ns = ete3.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1
        for node in t.traverse():
            node.set_style(ns)
            if node.is_leaf():
                if node.name.strip("'").title() in organisms:
                    node_face = ete3.TextFace(node.name.strip("'"), fsize=60, penwidth=10)
                    node.add_face(node_face, 1)
                else:
                    node_face = ete3.TextFace(node.name.strip("'"), fsize=45, penwidth=10, fgcolor="blue")
                    node.add_face(node_face, 1)

        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "radial",
                                   "{}.pdf".format(gene))
        t.render(destination, tree_style=ts)

    def _create_standard_gene_tree(self, newick_info, gene, organisms):
        """
        Creates the standard gene tree of each gene
        :param newick_info - newick-formatted tree representation
        :param gene - gene name
        :param organisms - list of organisms that contain the given gene
        """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.force_topology = True
        ts.title.add_face(ete3.TextFace(gene, fsize=120, bold=True), column=0)
        ts.show_leaf_name = False
        ts.show_scale = False
        ts.rotation = 90

        ns = ete3.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1
        for node in t.traverse():
            node.set_style(ns)
            if node.is_leaf():
                if node.name.strip("'").title() in organisms:
                    node_face = ete3.TextFace(node.name.strip("'"), fsize=60, penwidth=10)
                    node.add_face(node_face, 1)
                else:
                    node_face = ete3.TextFace(node.name.strip("'"), fsize=45, penwidth=10, fgcolor="blue")
                    node.add_face(node_face, 1)

        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "standard",
                                   "{}.pdf".format(gene))
        t.render(destination, tree_style=ts)


if __name__ == "__main__":
    mapper = Mapper()
    mapper.create_global_trees()
    mapper.create_gene_trees()
