import os

import matplotlib.pyplot as plt


class Mapper:
    """ Responsible for parsing the organisms taxon information and building plots of organisms classification """

    CLADE_IDX = 6

    def __init__(self):
        """ Constructor """
        self.organisms = self._read_organisms()

    def _read_organisms(self):
        """ Reads and constructs the organisms dictionary """
        result = {}
        path = os.path.join(os.getcwd(), "src", "data", "phylogeny", "taxa_info.txt")
        with open(path, "r") as f:
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
        txt_caption = "Taxonomic frequency of the analyzed dataset. The presented classes are:\n" \
                      "Mammalia - mammals; Actinopterygii - ray-finned fish; Archelosauria - turtles and archosaurs;" \
                      "\nLepidosauria - scaly reptiles; Amphibia - amphibians; Holocephali - cartilaginous fish; \n" \
                      "Coelacanthiformes - bony fish (ancient); Hyperoartia - jawless bony fish;\n" \
                      "Hyperotreti - hagfish; Phlebobranchia - sea squirts; Pterygota - winged insects; \n" \
                      "Rhabditina - nematodes."
        plt.figure(figsize=(8, 8))
        plt.xlabel(txt_caption)
        plt.ylim(0, max(y_values) + 5)
        plt.xticks(range(len(x_values)), x_values, rotation="45")
        plt.scatter(x_values, y_values)
        plt.savefig("taxa_info.pdf", quality=95, bbox_inches='tight')


if __name__ == "__main__":
    mapper = Mapper()
    mapper.classify()
