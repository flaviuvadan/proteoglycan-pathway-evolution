import matplotlib.pyplot as plt


class Mapper:
    """ Responsible for parsing the organisms taxon information and building plots of organisms classification """

    CLADE_IDX = 5

    def __init__(self):
        """ Constructor """
        self.organisms = self._read_organisms()

    def _read_organisms(self):
        """ Reads and constructs the organisms dictionary """
        result = {}
        with open("taxa_info.txt", "r") as f:
            f.readline()  # don't care about the title
            for line in f.readlines():
                split_line = line.split(",")
                clade = split_line[self.CLADE_IDX]
                if not result.get(clade):
                    result[clade] = 1
                else:
                    result[clade] = result.get(clade) + 1
        return result

    def classify(self):
        """ Creates a plot of the organisms classification """
        x_values = self.organisms.keys()
        y_values = self.organisms.values()
        plt.figure(figsize=(8, 10))
        fig_title = "Frequency of organisms classification based on bone, or cartilage, presence"
        plt.title(fig_title, fontsize=12)
        plt.xticks(range(len(x_values)), x_values, rotation="45")
        plt.yscale("log")
        plt.scatter(x_values, y_values)
        plt.savefig("taxa_info.pdf", quality=95)


if __name__ == "__main__":
    mapper = Mapper()
    mapper.classify()
