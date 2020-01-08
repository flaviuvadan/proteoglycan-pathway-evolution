import os

import matplotlib.pyplot as plt
import pandas as pd


class Distributor:
    """ Responsible for reading in the dnds data to create distributions """

    ORG_CLASS_IDX = 0
    ORG_NAME_IDX = 1

    # number of genes in this study
    NUM_GENES = 51

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
            final[sig_org] = sorted(histo_vals)
        return final

    def visualize(self):
        """ Creates a collection of distributions, one for each significant organism. The distributions represent the
        pattern in dnds ratios for each organism """
        fig, ax = plt.subplots(7, 3,
                               figsize=(12, 18))
        fig.suptitle("Significant organisms' dN/dS distributions",
                     fontsize=15)
        max_so_far = float('-inf')
        for idx, k in enumerate(self.histogram_data.keys()):
            plt.subplot(7, 3, idx + 1)
            fonts = {
                'fontsize': 11,
                'fontweight': 2,
                'verticalalignment': 'baseline',
                'horizontalalignment': 'center'
            }
            plt.title(" ".join(k.split("_")[:2]).capitalize(),
                      fontdict=fonts)
            data = self.histogram_data.get(k)
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
            max_so_far = max(max_so_far, max(self.histogram_data.get(k)))
            plt.hist(self.histogram_data.get(k), bins=21)
            plt.axvline(1,
                        color='k',
                        linestyle='dashed',
                        linewidth=1)
            plt.annotate("{} < 1\n{} >= 1".format(lt_1, gt_1), (5, 550))  # I hate arbitrary numbers like this, Starman
        title = os.path.join(os.getcwd(), "src", "data", "visualizations", "dnds", "histograms.pdf")
        print(max_so_far)
        plt.savefig(title,
                    format="pdf",
                    quality=95)


if __name__ == "__main__":
    dist = Distributor()
    dist.visualize()
