class Error(Exception):
    pass
class IllegalMassSizeError(Error):
    """Exception when size of mass vector is not 2^n"""
    def __init__(self, message):
        super(IlegalMassSizeError, self).__init__(message)
