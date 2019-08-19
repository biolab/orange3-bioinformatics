""" Bioinformatics add-on for Orange3 """

from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('orange3-bioinformatics').version
except DistributionNotFound:
    # package is not installed
    pass
