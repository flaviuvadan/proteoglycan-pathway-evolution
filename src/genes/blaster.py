import json
import os

from src import exceptions


class BlastKeys:
    OUTPUT = "BlastOutput2"
    REPORT = "report"
    RESULTS = "results"
    SEARCH = "search"
    QUERY_LENGTH = "query_len"
    HITS = "hits"
    HIT_SPECIES = "hsps"
    DESCRIPTION = "description"
    SCIENTIFIC_NAME = "sciname"
    TITLE = "title"
    E_VALUE = "evalue"
    QUERY_FROM = "query_from"
    QUERY_TO = "query_to"
    HIT_FROM = "hit_from"
    HIT_TO = "hit_to"
    IDENTITY = "identity"
    ALIGNMENT_LENGTH = "align_len"

    # percent query coverage for similarity filtering
    QUERY_COVER_LIMIT = 75


class Blaster:
    """ Responsible for parsing the gene trees BLAST results """
    SIG_DIGITS = 3
    ORG_TITLE_IDX = 0
    QRY_COV_IDX = 1
    E_VAL_IDX = 2
    IDENTITY_IDX = 3

    def __init__(self):
        self.blast_files = self._get_blast_file_paths()

    def _get_blast_file_paths(self):
        """ Returns a list of the paths to the BLAST files """
        path = os.path.join(os.getcwd(), "src", "data", "genes", "blast")
        files = os.listdir(path)
        paths = []
        for f in files:
            paths.append(os.path.join(path, f))
        return paths

    def _parse_blast_output(self, blast_output):
        """
        Parses the output of a given BLAST result
        :param dict blast_output: the BLAST result as a dictionary
        :return: tuple of results (org scientific name, % query cover, e-value, % identity)
        """
        output = blast_output.get(BlastKeys.OUTPUT)
        if not output:
            raise exceptions.BlastDictionaryAccessException("BLAST output not accessible")
        output = output[0]  # don't know why it's a list
        report = output.get(BlastKeys.REPORT)
        if not report:
            raise exceptions.BlastDictionaryAccessException("BLAST report not accessible")

        results = report.get(BlastKeys.RESULTS)
        if not results:
            raise exceptions.BlastDictionaryAccessException("BLAST results not accessible")

        search = results.get(BlastKeys.SEARCH)
        if not search:
            raise exceptions.BlastDictionaryAccessException("BLAST search results not accessible")

        hits = search.get(BlastKeys.HITS)
        if not hits:
            raise exceptions.BlastDictionaryAccessException("BLAST hits not accessible")

        results = []
        query_len = search.get(BlastKeys.QUERY_LENGTH)
        # there can be multiple hits if BLAST decides to use the complete genome sequence of an organism
        # for multiple hits, we only care about the scores of the ones with high query coverage
        # so we set the return values according to that
        for hit in hits:
            query_coverage = 0
            identity = 0
            align_total = 0
            e_val = float('-inf')
            hsps = hit.get(BlastKeys.HIT_SPECIES)
            for hsp in hsps:
                curr_query_diff = hsp.get(BlastKeys.QUERY_TO) - hsp.get(BlastKeys.QUERY_FROM)
                curr_query = round(curr_query_diff / query_len, self.SIG_DIGITS) * 100
                # we set the e-value, identity, and title according to query coverage
                if curr_query > query_coverage:
                    query_coverage = curr_query
                    curr_identity = round(hsp.get(BlastKeys.IDENTITY) / hsp.get(BlastKeys.ALIGNMENT_LENGTH),
                                          self.SIG_DIGITS) * 100
                    identity = max(identity, curr_identity)
                    align_total += hsp.get(BlastKeys.ALIGNMENT_LENGTH)
                    e_val = max(e_val, hsp.get(BlastKeys.E_VALUE))

            # if we're here, it's fine to process the dict, it'll have something, no exception
            hit_description = hit.get(BlastKeys.DESCRIPTION)[0]  # not sure why this is a list
            org_title = hit_description.get(BlastKeys.TITLE)
            if query_coverage < BlastKeys.QUERY_COVER_LIMIT:
                continue
            results.append((org_title, query_coverage, e_val, identity))
        return results

    def parse_blast_files(self):
        """ Parses the BLAST files and builds """
        blast_results_file = os.path.join(os.getcwd(), "src", "data", "genes", "blast_results.txt")
        with open(blast_results_file, "w") as brf:
            # make a nice header, otherwise it's hard to understand the format b/c of long titles
            brf.write("GENE NAME\nORGANISM TITLE\nQUERY COVERAGE (%), E-VALUE, IDENTITY (%)\n\n")

            for blast_file in self.blast_files:
                gene_name = blast_file.split("/")[-1].split(".")[0]
                with open(blast_file, "r") as bf:
                    blast = json.load(bf)
                    try:
                        orgs = self._parse_blast_output(blast)
                    except exceptions.BlastDictionaryAccessException:
                        continue  # can safely skip
                    if not orgs:
                        continue
                    else:
                        brf.write("{}\n".format(gene_name))
                        for o in orgs:
                            brf.write("{}\n".format(o[self.ORG_TITLE_IDX]))
                            vals = "{}%, {}, {}%\n".format(o[self.QRY_COV_IDX], o[self.E_VAL_IDX], o[self.IDENTITY_IDX])
                            brf.write(vals)
                        brf.write("\n")


if __name__ == "__main__":
    blaster = Blaster()
    blaster.parse_blast_files()
