""" Utility functions for GEO module """
import os
import io

from re import compile
from numpy import isnan, median, array, asarray
from tempfile import NamedTemporaryFile
from urllib.request import urlopen
from Orange.data import Table

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.geo.config import DOMAIN, FTP_DIR, FTP_NCBI

unknown = float('nan')
attribute_reg = compile("![A-Za-z]*_([A-Za-z_]*) = (.*)$")


def create_table(domain, X, Y, metas):
    class_var = domain.class_var
    meta_atts = domain.metas

    if Y:
        Y = array([[class_var.to_val(row)] for row in Y], dtype=float)
    if metas:
        metas = array([[c.to_val(v) for c, v in zip(meta_atts, row)] for row in metas], dtype=object)

    return Table(domain, asarray(X), Y=Y, metas=metas)


def spots_mean(x):
    vs = [v for v in x if not isnan(v)]
    if len(vs) == 0:
        return unknown
    else:
        return sum(vs) / len(vs)


def spots_median(x):
    vs = [v for v in x if not isnan(v)]
    if len(vs) == 0:
        return unknown
    return median(vs)


def spots_min(x):
    vs = [v for v in x if not isnan(v)]
    if len(vs) == 0:
        return unknown
    else:
        return min(vs)


def spots_max(x):
    vs = [v for v in x if not isnan(v)]
    if len(vs) == 0:
        return unknown
    else:
        return max(vs)


def to_float(value):
    return float(value) if value != 'null' else unknown


def parse_attribute_line(string):
    return attribute_reg.search(string).groups()


def gds_is_cached(gds_name):
    return os.path.isfile(
        os.path.join(serverfiles.localpath(DOMAIN), gds_name + ".soft.gz"))


def sniff_size(file_obj):
    if isinstance(file_obj, io.FileIO):
        return os.fstat(file_obj.fileno()).st_size
    return None


def gds_download_url(gds_name):
    """ Return the download url for a GDS id `gdsname`.
    """
    return "ftp://{}/{}/{}.soft.gz".format(FTP_NCBI, FTP_DIR, gds_name)


def gds_ensure_downloaded(gds_name, progress=None):
    """ Ensure the GDS dataset is available locally in cache.
    """
    if gds_is_cached(gds_name):
        return
    else:
        gds_download(gds_name, progress=progress)


def gds_download(gds_name, progress=None):
    """ Download the GDS dataset into the cache.
    """
    gds_url = gds_download_url(gds_name)
    basename = gds_name + ".soft.gz"
    target_path = os.path.join(serverfiles.localpath(DOMAIN), basename)

    temp = NamedTemporaryFile(prefix=basename + "-", dir=serverfiles.localpath(DOMAIN), delete=False)
    try:
        retrieve_url(gds_url, temp, progress=progress)
    except BaseException as err:
        try:
            temp.close()
            os.remove(temp.name)
        except (OSError, IOError):
            pass
        raise err
    else:
        temp.close()
        os.replace(temp.name, target_path)


def retrieve_url(url, dest_obj, progress=None):
    """ Retrieve an `url` writing it to an open file-like `destobj`.

    Args:
        url: The source url.
        dest_obj: An file-like object opened for writing.

        progress: (int, int) -> None optional
                  An optional progress callback function. Will be called
                  periodically with `(transfered, total)` bytes count. `total`
                  can be `-1` if the total contents size cannot be
                  determined beforehand.

    Returns:

    """

    with urlopen(url, timeout=10) as stream:
        length = stream.headers.get("content-length", None)
        if length is not None:
            length = int(length)
        copyfileobj(stream, dest_obj, size=length, progress=progress)


def copyfileobj(src, dst, buffer=2 ** 15, size=None, progress=None):
    """ Like shutil.copyfileobj but with progress reporting.

    Args:
        src: Source file object
        dst: Destination file object
        buffer: Buffer size
        size: Total `src` contents size if available.

        progress: : (int, int) -> None optional
                    An optional progress callback function. Will be called
                    periodically with `(transferred, total)` bytes count. `total`
                    can be `-1` if the total contents size cannot be
                    determined beforehand.

    Returns:

    """

    count = 0
    if size is None:
        size = sniff_size(src)

    while True:
        data = src.read(buffer)
        dst.write(data)
        count += len(data)
        if progress is not None:
            progress(count, size if size is not None else -1)
        if not data:
            break

    if size is None and progress is not None:
        progress(count, count)

    return count
