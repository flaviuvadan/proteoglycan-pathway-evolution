import os

from ete3 import Tree, TreeStyle, AttrFace


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
        """ Creates the radial version of the global tree """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = 'c'
        ts.arc_span = 360
        ts.force_topology = True
        ts.legend = None
        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "radial", "global.pdf")
        t.render(destination, tree_style=ts)

    def _create_standard_tree(self, newick_info):
        """ Creates the standard version of the global tree """
        t = Tree(newick_info, format=1)  # see ete3/coretype/TreeNode doc for format
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.force_topology = True
        ts.legend = None
        ts.rotation = 90
        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "standard", "global.pdf")
        t.render(destination, tree_style=ts)

    def create_gene_trees(self):
        """ Creates the 51 trees, one per gene """
        with open(os.path.join(os.getcwd(), "src", "data", "phylogeny", "phyliptree.phy"), "r") as f:
            newick_info = f.read()

    def _create_radial_gene_tree(self, newick_info):
        """ Creates the radial gene tree of each gene """


if __name__ == "__main__":
    mapper = Mapper()
    mapper.create_global_trees()
