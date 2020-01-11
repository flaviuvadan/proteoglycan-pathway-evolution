class OrganismBoneGroups:
    """ A namespace class that holds how organisms should be grouped based on bone presence """
    BONE_CART = "bone_cartilage"
    CART_ONLY = "cartilage"
    NO_BONE_NO_CART = "neither"

    def __init__(self):
        """ Constructor """
        self.org_class_map = self._get_org_class_map()
        self.class_org_map = self._get_class_org_map()

    def _get_org_class_map(self):
        """ Creates a mapping of organism genus and species to its bone class """
        return {
            "homo_sapiens": self.BONE_CART,
            "mus_musculus": self.BONE_CART,
            "ornithorhynchus_anatinus": self.BONE_CART,
            "danio_rerio": self.BONE_CART,
            "oryzias_latipes_hsok": self.CART_ONLY,
            "takifugu_rubripes": self.CART_ONLY,
            "lepisosteus_oculatus": self.CART_ONLY,
            "anas_platyrhynchos_platyrhynchos": self.BONE_CART,
            "gallus_gallus": self.BONE_CART,
            "chrysemys_picta_bellii": self.BONE_CART,
            "crocodylus_porosus": self.BONE_CART,
            "notechis_scutatus": self.BONE_CART,
            "anolis_carolinensis": self.BONE_CART,
            "xenopus_tropicalis": self.BONE_CART,
            "callorhinchus_milii": self.CART_ONLY,
            "latimeria_chalumnae": self.BONE_CART,
            "petromyzon_marinus": self.BONE_CART,
            "eptatretus_burgeri": self.BONE_CART,
            "ciona_intestinalis": self.NO_BONE_NO_CART,
            "drosophila_melanogaster": self.NO_BONE_NO_CART,
            "caenorhabditis_elegans": self.NO_BONE_NO_CART,
        }

    def _get_class_org_map(self):
        """ Creates a mapping of bone class to organisms """
        final = {}
        for k, v in self.org_class_map.items():
            if not final.get(v):
                final[v] = [k]
            else:
                final[v].append(k)
        return final


class OrganismHabitatGroups:
    """ A namespace class that holds how organisms should be grouped based on habitat """

    TERR = "terrestrial"
    AQUA = "aquatic"
    TERR_AQUA = "terrestrial/aquatic"

    def __init__(self):
        """ Constructor """
        self.org_class_map = self._get_org_class_map()
        self.class_org_map = self._get_class_org_map()

    def _get_org_class_map(self):
        """ Creates a mapping of organism genus and species to its habitat class """
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

    def _get_class_org_map(self):
        """ Creates a mapping of org class to species """
        final = {}
        for k, v in self.org_class_map.items():
            if not final.get(v):
                final[v] = [k]
            else:
                final[v].append(k)
        return final
