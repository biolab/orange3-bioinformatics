""" Bioinformatics add-on for Orange3 """

from importlib.metadata import distribution, PackageNotFoundError

try:
    __version__ = distribution('orange3-bioinformatics').version
except PackageNotFoundError:
    # package is not installed
    pass
