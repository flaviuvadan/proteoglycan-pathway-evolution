import os

import ete3 as ete


class Mapper:
    """ Mapper is reponsible for creating phylogenetic trees for all the significant organisms that were selected to be
    the subset of all the organisms whose phylogenetic relationship have been created by mapper for all the organisms"""

    def create_tree(self):
        """ Creates the main tree """
        path = os.path.join(os.getcwd(), "src", "data", "phylogeny", "sig_orgs_phyliptree.phy")
        with open(path, "r") as f:
            newick_info = f.read()
            self._create_tree(newick_info)

    def _create_tree(self, newick_info):
        """
        Creates the radial tree that represents the phylogeny as dictated by the given info
        :param newick_info: newick-formatted tree representation
        """
        t = ete.Tree(newick_info, format=1)  # see ete/coretype/TreeNode doc for format
        ts = ete.TreeStyle()
        ts.mode = 'c'
        ts.arc_span = 360
        ts.show_leaf_name = False
        ts.force_topology = True
        ts.show_scale = False
        ts.scale = 120
        tree_face = ete.TextFace("Significant organisms phylogenetic relationship", fsize=100, bold=True)
        ts.title.add_face(tree_face, column=0)

        ns = ete.NodeStyle()
        ns["hz_line_type"], ns["vt_line_type"] = 0, 0
        ns["hz_line_color"], ns["vt_line_color"] = "black", "black"
        ns["hz_line_width"], ns["vt_line_width"] = 1, 1

        for node in t.traverse():
            node.set_style(ns)
            if not node.is_leaf():
                face = ete.TextFace(node.name, fsize=40, penwidth=10, fgcolor="blue")
                node.add_face(face, 1)
            else:
                face = ete.TextFace(node.name, fsize=65, penwidth=10)
                node.add_face(face, 1)
        destination = os.path.join(os.getcwd(), "src", "data", "visualizations", "phylogeny", "sig_orgs_phylo.pdf")
        t.render(destination, tree_style=ts, dpi=500)


if __name__ == "__main__":
    mapper = Mapper()
    mapper.create_tree()
