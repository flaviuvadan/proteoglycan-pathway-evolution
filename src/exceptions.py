class PGException(Exception):
    """ A general exception for the proteoglycan pathway project """
    pass


class OrthologRequestException(PGException):
    """ An exception used for handling Ensembl API non-200 responses """
    pass


class EmptyOrthologData(PGException):
    """ An exception used for handling the emptiness of the data field of the Ensembl response """
    pass


class EmptyHomologyInformation(PGException):
    """ An exception used for handling empty homology information - this should rarely be the case, if ever """
    pass


class EmptySquence(PGException):
    """ An exception used for handling empty sequences in the Ensembl API response """
    pass


class SpeciesRequestException(PGException):
    """ An exception used for handling GBIF API non-200 responses """
    pass


class ProjectedBuildOrganismException(PGException):
    """ An exception used for indicating an organism has a projected built genome in Ensembl """
    pass


class BlastDictionaryAccessException(PGException):
    """ An exception used for indicating that a BLAST dictionary does not have something stored at a passed key """
    pass


class EmptyNamespaceClassException(PGException):
    """ An exception used for indicating that no namespace class has been passed for binning organisms """
    pass
