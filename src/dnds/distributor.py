import os

import matplotlib.pyplot as plt
import pandas as pd
import statistics


class OrganismGroups:
    """ A namespace class that holds how organisms should be grouped """

    TERR = "terrestrial"
    AQUA = "aquatic"
    TERR_AQUA = "terrestrial/aquatic"

    def __init__(self):
        """ Constructor """
        self.org_class_map = self._get_org_class_mapping()
        self.class_org_map = self._get_class_org_mapping()

    def _get_org_class_mapping(self):
        """ Creates a mapping of organism genus and species to its class """
        return {
            "homo_sapiens": self.TERR,
            "mus_musculus": self.TERR,
            "ornithorhynchus_anatinus": self.TERR_AQUA,
            "danio_rerio": self.AQUA,
            "oryzias_latipes_hsok": self.AQUA,
            "takifugu_rubripes": self.AQUA,
            "lepisosteus_oculatus": self.AQUA,
            "anas_platyrhynchos_platyrhynchos": self.TERR_AQUA,
            "gallus_gallus": self.TERR,
            "chrysemys_picta_bellii": self.AQUA,
            "crocodylus_porosus": self.TERR_AQUA,
            "notechis_scutatus": self.TERR,
            "anolis_carolinensis": self.TERR,
            "xenopus_tropicalis": self.TERR_AQUA,
            "callorhinchus_milii": self.AQUA,
            "latimeria_chalumnae": self.AQUA,
            "petromyzon_marinus": self.AQUA,
            "eptatretus_burgeri": self.AQUA,
            "ciona_intestinalis": self.AQUA,
            "drosophila_melanogaster": self.TERR,
            "caenorhabditis_elegans": self.TERR_AQUA,
        }

    def _get_class_org_mapping(self):
        """ Creates a mapping of org class to species """
        final = {}
        for k, v in self.org_class_map.items():
            if not final.get(v):
                final[v] = [k]
            else:
                final[v].append(k)
        return final


class Distributor:
    """ Responsible for reading in the dnds data to create distributions """

    ORG_CLASS_IDX = 0
    ORG_NAME_IDX = 1

    # number of genes in this study
    NUM_GENES = 51

    def __init__(self):
        """ Constructor """
        self.sig_orgs = self._get_significant_orgs()
        self.dnds_data = self._get_hg_data()

    def _get_significant_orgs(self):
        """ Reads in the significant organisms selected for this study """
        orgs = {}
        path = os.path.join(os.getcwd(), "src", "data", "phylogeny", "significant_organisms.txt")
        with open(path, "r") as f:
            f.readline()  # don't care about the top line
            for line in f.readlines():
                org_name = line.split(",")[self.ORG_NAME_IDX]
                org = "_".join(org_name.lower().split())
                orgs[org] = 1
        return orgs

    def _get_hg_data(self):
        """ Creates a dictionary of the data to be plotted """
        final = {}
        for sig_org in self.sig_orgs.keys():
            sig_org_csv = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(sig_org))
            df = pd.read_csv(sig_org_csv)
            histo_vals = []
            for x in df.values:
                histo_vals = histo_vals + list(x[1:])
            final[sig_org] = sorted(histo_vals)
        return final

    def visualize(self):
        """ Creates a histogram of the independent significant organisms dN/dS distribution along with a histogram of
        the binned version for visualizing the distribution per group """
        self._visualize_independent_orgs()
        self._visualize_grouped_orgs()

    def _visualize_independent_orgs(self):
        """ Creates a collection of distributions, one for each significant organism. The distributions represent the
        pattern in dnds ratios for each organism """
        fig, ax = plt.subplots(7, 3,
                               figsize=(12, 18))
        fig.suptitle("Significant organisms' dN/dS distributions",
                     fontsize=15)
        for idx, k in enumerate(self.dnds_data.keys()):
            plt.subplot(7, 3, idx + 1)
            fonts = {
                'fontsize': 11,
                'fontweight': 2,
                'verticalalignment': 'baseline',
                'horizontalalignment': 'center'
            }
            plt.title(" ".join(k.split("_")[:2]).capitalize(),
                      fontdict=fonts)
            data = self.dnds_data.get(k)
            lt_1 = len([x for x in data if x < 1])
            gt_1 = len(data) - lt_1
            plt.axis([0, 7, 0, 800])
            # have to remove indices manually since I cannot figure out why sharex and sharey have no effect when used
            # with plt.subplots
            no_x_tick = [x for x in range(1, 19)]
            if idx + 1 in no_x_tick:
                plt.xticks([])
            y_tick = [y for y in range(1, 20, 3)]
            if idx + 1 not in y_tick:
                plt.yticks([])
            plt.hist(self.dnds_data.get(k), bins=21)
            plt.axvline(1,
                        color='k',
                        linestyle='dashed',
                        linewidth=1)
            # I hate arbitrary numbers like this, Starman
            plt.annotate("{} < 1\n{} >= 1".format(lt_1, gt_1), (5, 550))
        title = os.path.join(os.getcwd(), "src", "data", "visualizations", "dnds", "grouped_orgs", "histograms.pdf")
        plt.savefig(title,
                    format="pdf",
                    quality=95)

    def _visualize_grouped_orgs(self):
        """ Creates a collection of distributions for each organism bin. The bins are: bone, cartilage, neither,
        terrestrial, aquatic, both """
        binned = {}
        groups = OrganismGroups()
        x_labels = None  # this should only be set once, we need to preserve the order of the genes
        for org_class, orgs in groups.class_org_map.items():
            binned[org_class] = {}
            for org in orgs:
                org_csv = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(org))
                df = pd.read_csv(org_csv)
                if not x_labels:
                    # first one is "organism", should be set on the first try
                    x_labels = list(df.columns[1:])
                for gene in x_labels:  # after all, those are genes...
                    if not binned[org_class].get(gene):
                        binned[org_class][gene] = list(df[gene])
                    else:
                        binned[org_class][gene].extend(list(df[gene]))
        # go over each gene, grab average for each, plot, move on
        # subplots 111
        # +-=0.2 width for plots
        x_range = list(range(1, len(x_labels) + 1))
        plt.figure(figsize=(10, 10))
        interval_width = 0.5
        for idx, gene in enumerate(x_labels):
            plt.bar(list(map(lambda x: x - interval_width, x_range)), statistics.mean(binned[groups.TERR][gene]),
                    width=interval_width,
                    color='b',
                    align='center')
            plt.bar(x_range, statistics.mean(binned[groups.AQUA][gene]),
                    width=interval_width,
                    color='g',
                    align='center')
            plt.bar(list(map(lambda x: x + interval_width, x_range)), statistics.mean(binned[groups.TERR_AQUA][gene]),
                    width=interval_width,
                    color='r',
                    align='center')
        plt.xticks(range(1, len(x_labels) + 1), x_labels, rotation=90)
        plt.show()


if __name__ == "__main__":
    dist = Distributor()
    dist._visualize_grouped_orgs()
