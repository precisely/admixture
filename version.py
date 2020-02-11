# -*- coding: utf-8 -*-
"""
This module holds the version information so it only has to be changed
in one place. Based on `<http://bit.ly/16LbuJF>`_
"""


def _safe_int(string):
    """ Simple function to convert strings into ints without dying. """
    try:
        return int(string)
    except ValueError:
        return string


__version__ = '0.0.0'
VERSION = tuple(_safe_int(x) for x in __version__.split('.'))

