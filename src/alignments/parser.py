import os

from matplotlib import collections as matcoll
import matplotlib.pyplot as plt

from src.alignments import constants

plt.style.use('seaborn-white')


class Parser:
    """ Parser is responsible for parsing the multiple sequence alignment files from data/alignments """

    GENE_NAME_IDX = 0  # index of the gene name when parsing gene_filename.txt
    ALL = "all_genes"  # whether to plot all genes when using build_frequency_plots_by_function

    def __init__(self):
        """ Constructor """
        self.alignment_file_paths = self._build_alignment_paths()
        self.alignments = self._parse_alignments()
        self.frequencies = self._build_frequencies()
        self.frequencies_by_function = self._build_frequencies_by_function()

    def _get_alignment_filenames(self):
        """ Returns a list of all the MSA filenames """
        return os.listdir(os.path.join(os.getcwd(), "src", "data", "alignments"))

    def _build_alignment_paths(self):
        """ Returns a list of all the alignments full file paths """
        alignment_file_names = self._get_alignment_filenames()
        return [os.path.join(os.getcwd(), "src", "data", "alignments", afn) for afn in alignment_file_names]

    def _parse_alignments(self):
        """ Parses the MSAs to clean up the 60 character break added by ClustalW """
        final = {}
        for afp in self.alignment_file_paths:
            filename = afp.split("/")[-1]  # the filename.txt is sufficient
            if not final.get(filename):
                final[filename] = {}
            with open(afp, "r") as f:
                org = ""  # helps keep track of the last organism
                for line in f.readlines():
                    if not line.strip():
                        continue  # it can happen that we have empty lines
                    line = line.strip("\n")
                    if ">" in line:
                        org = line.strip(">")
                        if not final.get(filename).get(org):
                            final[filename][org] = ""
                    else:
                        final[filename][org] = final.get(filename).get(org) + line
        return final

    def _build_frequencies(self):
        """ Builds a list of lists containing tuples indicating the number of characters per alignment column. This is
        used for building a simple visualization of the generated MSAs """
        freqs = {}
        for gene in self.alignments.keys():
            if not freqs.get(gene):
                freqs[gene] = []
            for i in range(len(list(self.alignments.get(gene).values())[0])):  # all alignments have the same length
                freq_i = 0
                for seq in self.alignments.get(gene).values():
                    if seq[i] != "-":
                        freq_i = freq_i + 1
                freqs[gene].append((i, freq_i))
        return freqs

    def _build_frequencies_by_function(self):
        """ Builds a mapping of gene functions to gene frequencies used for grouping visualizations in MSAs """
        groups = {}
        for gene in self.frequencies.keys():
            g = gene.split("_")[self.GENE_NAME_IDX]
            group = constants.GENE_FUNCTIONS.get(g)
            if not groups.get(group):
                groups[group] = {}
            groups[group][gene] = self.frequencies.get(gene)
        return groups

    def build_block_alignments(self):
        """ Reads the parsed multiple sequence alignments and prints to a file in block format """
        for align_key in self.alignments.keys():
            # I do not remember why I've put this in src/data/organisms... -.-
            path = os.path.join(os.getcwd(), "src", "data", "organisms", align_key)
            with open(path, "w") as f:
                for org_key in self.alignments.get(align_key).keys():
                    f.write(">{}\n".format(org_key))
                    f.write("{}\n".format(self.alignments.get(align_key).get(org_key).replace("-", "")))

    def build_frequency_plots(self, subset=False, subset_num=0):
        """
        Builds the frequency plots associated with the MSAs
        :param subset: whether to respect the subset_num parameter
        :param subset_num: the number of genes to plot as a subset of the total
        """
        num_plots = len(list(self.frequencies.keys()))
        fig = plt.figure(figsize=(8, 40))  # these were figured out by trial and error
        fig.subplots_adjust(hspace=2)
        for idx, gene in enumerate(list(self.frequencies.keys())):
            if subset and idx == subset_num:
                break
            values = self.frequencies.get(gene)
            x_vals = [v[0] for v in values]
            y_vals = [v[1] for v in values]
            lines = []
            for i in range(len(x_vals)):
                pair = [(x_vals[i], 0), (x_vals[i], y_vals[i])]
                lines.append(pair)
            line_col = matcoll.LineCollection(lines)
            ax = fig.add_subplot(num_plots, 1, idx + 1)
            ax.add_collection(line_col)
            plt.scatter(x_vals, y_vals, marker="None")
            plt.xlim(0, len(lines))
            plt.yticks([])
            title = "{}".format(gene.split('_')[self.GENE_NAME_IDX])
            plt.title(title)
        path = "src/data/visualizations/MSAs.pdf"
        plt.savefig(path, format="pdf", quality=95, bbox_inches='tight')

    def build_frequency_plots_by_function(self):
        """ Builds the frequency plots associated with the MSAs. Plots are grouped
        based on gene function """
        for i, func in enumerate(self.frequencies_by_function):
            num_plots = len(list(self.frequencies_by_function.get(func)))
            fig = plt.figure(figsize=(8, num_plots * .8))
            fig.subplots_adjust(hspace=2)
            for idx, gene in enumerate(list(self.frequencies_by_function.get(func).keys())):
                values = self.frequencies_by_function.get(func).get(gene)
                x_vals = [v[0] for v in values]
                y_vals = [v[1] for v in values]
                lines = []
                for i in range(len(x_vals)):
                    pair = [(x_vals[i], 0), (x_vals[i], y_vals[i])]
                    lines.append(pair)
                line_col = matcoll.LineCollection(lines)
                ax = fig.add_subplot(num_plots, 1, idx + 1)
                ax.add_collection(line_col)
                plt.scatter(x_vals, y_vals, marker="None")
                plt.xlim(0, len(lines))
                plt.yticks([])
                title = "{}".format(gene.split('_')[self.GENE_NAME_IDX])
                plt.title(title)
            path = "src/data/visualizations/alignments/{}_msa.pdf".format(func.lower())
            plt.savefig(path, format="pdf", quality=95, bbox_inches='tight')


if __name__ == "__main__":
    parser = Parser()
    # parser.build_block_alignments()
    # parser.build_frequency_plots()
    parser.build_frequency_plots_by_function()
