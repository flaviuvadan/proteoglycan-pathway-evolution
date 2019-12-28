class Mapper:
    """ Responsible for finding similarities between phylogenetic trees and creating visualizations for the results """

    def __init__(self):
        self.tree_data = self._load_tree_data()
        self.similarities = self._find_similarities()

    def _find_similarities(self):
        return 1

    def _load_tree_data(self):
        _ = self._get_path_to_tree_data()
        return 1

    def _get_path_to_tree_data(self):
        return 1

    def create_similarity_trees(self):
        pass


if __name__ == "__main__":
    mapper = Mapper()
    mapper.create_similarity_trees()
