"""

"""
from __future__ import absolute_import

from . import conf

REST_API = "http://rest.kegg.jp/"


def slumber_service():
    """
    Return a rest based service using `slumber` package
    """
    import slumber

    if not hasattr(slumber_service, "_cached"):

        class DecodeSerializer(slumber.serialize.BaseSerializer):
            key = "decode"
            content_types = ["text/plain"]

            def loads(self, data):
                return data

        # for python 2/3 compatibility
        serializer = slumber.serialize.Serializer(default="decode", serializers=[DecodeSerializer()])
        slumber_service._cached = slumber.API(REST_API, serializer=serializer)
    return slumber_service._cached


default_service = slumber_service

web_service = slumber_service
