import os

import matplotlib.pyplot as plt
import pandas as pd


class Distributor:
    """ Responsible for reading in the dnds data to create distributions """

    ORG_CLASS_IDX = 0
    ORG_NAME_IDX = 1

    def __init__(self):
        """ Constructor """
        self.sig_orgs = self._get_significant_orgs()
        self.histogram_data = self._get_hg_data()

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
            final[sig_org] = histo_vals
        return final

    def visualize(self):
        """ Creates a collection of distributions, one for each significant organism. The distributions represent the
        pattern in dnds ratios for each organism """
        fig, ax = plt.subplots(4, 3, figsize=(10, 10))
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle("Significant organisms' dN/dS distributions",
                     fontsize=15)
        for idx, k in enumerate(self.histogram_data.keys()):
            plt.subplot(4, 3, idx + 1)
            fonts = {
                'fontsize': 11,
                'fontweight': 2,
                'verticalalignment': 'baseline',
                'horizontalalignment': 'center'
            }
            plt.title(" ".join(k.split("_")).capitalize(), fontdict=fonts)
            plt.hist(self.histogram_data.get(k), bins=30)
        title = os.path.join(os.getcwd(), "src", "data", "visualizations", "dnds", "histograms.pdf")
        plt.savefig(title,
                    format="pdf",
                    quality=95)


if __name__ == "__main__":
    dist = Distributor()
    dist.visualize()
