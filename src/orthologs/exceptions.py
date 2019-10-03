class OrthologException(Exception):
    """ A general exception for the orthologs package """
    pass


class OrthologRequestException(OrthologException):
    """ An exception used for handling Ensembl API non-200 calls """
