import os


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
        # final result structure should be:
        #   alignments = {filename: {alignments}}
        final = {}
        for afp in self.alignment_file_paths:
            filename = afp.split("/")[-1]  # the filename.txt is sufficient
            if not final.get(filename):
                final[filename] = {}
            with open(afp, "r") as f:
                for line in f.readlines():
                    split_line = line.split()
                    if not split_line:  # this can happen because we have line breaks between alignment blocks
                        continue
                    org, alignment = split_line[self.ORGANISM_IDX], split_line[self.ALIGNMENT_IDX]
                    if not final.get(filename).get(org):
                        final[filename][org] = ""
                    final[filename][org] = final.get(filename).get(org) + alignment
        return final

    def _build_frequencies(self):
        """ Builds a list of lists containing tuples indicating the number of characters per alignment column. This is
        used for building a simple visualization of the generated MSAs """
        freqs = {}
        for gene in self.alignments.keys():
            if not freqs.get(gene):
                freqs[gene] = {}
            gene_alignment = self.alignments.get(gene)



if __name__ == "__main__":
    parser = Parser(os.getcwd())
