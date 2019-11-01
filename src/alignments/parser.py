import os

from matplotlib import collections as matcoll
import matplotlib.pyplot as plt

plt.style.use('seaborn-white')


class Parser:
    """ Parser is responsible for parsing the multiple sequence alignment files from data/alignments """

    ORGANISM_IDX = 0
    ALIGNMENT_IDX = 1

    def __init__(self, base_path):
        """ Constructor """
        self.base_path = base_path
        self.alignment_file_paths = self._build_alignment_paths()
        self.alignments = self._parse_alignments()
        self.frequencies = self._build_frequencies()

    def _get_alignment_filenames(self):
        """ Returns a list of all the MSA filenames """
        return os.listdir(self.base_path)

    def _build_alignment_paths(self):
        """ Returns a list of all the alignments full file paths """
        alignment_file_names = self._get_alignment_filenames()
        return [os.path.join(self.base_path, afn) for afn in alignment_file_names]

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
            for i in range(len(list(self.alignments.get(gene).values())[0])):
                freq_i = 0
                for seq in self.alignments.get(gene).values():
                    try:
                        if seq[i] != "-":
                            freq_i = freq_i + 1
                    except Exception as e:
                        pass
                freqs[gene].append((i, freq_i))
        return freqs

    def build_frequency_plots(self):
        """ Builds the frequency plots associated with the MSAs """
        num_plots = len(list(self.frequencies.keys()))
        fig = plt.figure(figsize=(10, 50))
        fig.subplots_adjust(hspace=1)
        for idx, gene in enumerate(list(self.frequencies.keys())):
            values = self.frequencies.get(gene)
            x_vals = [x[0] for x in values]
            y_vals = [x[1] for x in values]
            lines = []
            for i in range(len(x_vals)):
                pair = [(x_vals[i], 0), (x_vals[i], y_vals[i])]
                lines.append(pair)
            linecoll = matcoll.LineCollection(lines)
            ax = fig.add_subplot(num_plots, 1, idx + 1)
            ax.add_collection(linecoll)
            plt.scatter(x_vals, y_vals, marker="None")
            plt.yticks([])
            title = "{} - MSA".format(gene.split("_")[0])
            plt.title(title)
        plt.savefig("MSAs.pdf", format="pdf", quality=95)


if __name__ == "__main__":
    parser = Parser(os.getcwd())
    parser.build_frequency_plots()
