""" minimal environ module for serverfiles to work """

import os
import sys


def _settings_dir():
    name = "Orange3"
    if sys.platform == "darwin":
        return os.path.join(
            os.path.expanduser("~/Library/Application Support/"), name)
    elif sys.platform == "win32":
        return os.path.join(
            os.environ.get("APPDATA", os.path.expanduser("~/")), name)
    else:
        # TODO: use XGD_* env vars when available
        return os.path.expanduser("~/." + name)


def _cache_dir():
    name = "Orange3"
    if sys.platform == "darwin":
        return os.path.join(os.path.expanduser("~/Library/Caches"), name)
    else:
        settings_dir = _settings_dir()
        return os.path.join(settings_dir, "buffer")


buffer_dir = _cache_dir()
