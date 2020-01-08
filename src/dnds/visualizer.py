import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn

from src.dnds import loader


class Visualizer(loader.Loader):
    """ Responsible for creating dn/ds visualizations for all the studied genes """

    def __init__(self):
        """ Constructor """
        super(Visualizer, self).__init__()
        self._create_heatmap_csv()

    def _create_heatmap_csv(self):
        """ Creates a heatmap CSV file for visualizing the results """
        for sig_org in self.scores.keys():
            sig_org_file_path = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(sig_org))
            with open(sig_org_file_path, "w") as sof:
                header = "organism,"
                # each org will have itself in the dnds dictionary
                # we only want this to build the header once
                for k in self.scores[sig_org][sig_org].keys():
                    header += "{},".format(k)
                header = header[:len(header) - 1]  # remove the last ","
                header += "\n"
                sof.write(header)
                for org in self.scores[sig_org].keys():
                    if sig_org == org:  # no point in adding the same one
                        continue
                    vals = ",".join([str(x) for x in self.scores[sig_org][org].values()])
                    line = "{},{}\n".format(org, vals)
                    sof.write(line)

    def visualize(self):
        """ Creates the 2D dn/ds plots for each gene """
        plt.figure(figsize=(20, 20))
        for sig_org in self.scores.keys():
            sig_org_csv = os.path.join(os.getcwd(), "src", "data", "dnds", "{}.csv".format(sig_org))
            df = pd.read_csv(sig_org_csv)
            x_vals = df.columns[1:]
            y_vals = df.organism
            y_vals_clean = [" ".join(y.split("_")).capitalize() for y in y_vals]
            df = df.drop(columns="organism")
            plt.title("{} dN/dS ratios".format(" ".join(sig_org.split("_")).capitalize()), fontsize=30)
            sn.heatmap(df,
                       vmin=0,
                       vmax=5,
                       cbar=True,
                       xticklabels=x_vals,
                       yticklabels=y_vals_clean,
                       cmap="Blues",
                       linewidths=0.5)
            save_fig_path = os.path.join(os.getcwd(), "src", "data", "visualizations", "dnds", "{}.pdf".format(sig_org))
            plt.savefig(save_fig_path,
                        quality=95,
                        bbox_inches='tight')
            plt.clf()


if __name__ == "__main__":
    visualizer = Visualizer()
    visualizer.visualize()
