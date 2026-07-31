#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# File   : __init__.py
# Author : Honghua Dong
# Email  : dhh19951@gmail.com
#
# Distributed under terms of the MIT license.

from .globalpool import *

# GOBoost's published configurations use max_height=1 and therefore only need
# global pooling.  The optional hierarchical poolers pull in PyG extension
# modules that are not used for inference and make a portable LAFA image much
# harder to build.  Keep them available when their optional dependencies are
# installed, without making them mandatory for ordinary imports.
try:
    from .edgepool_dgl import EdgePooling
    from .louvain_dgl import LouvainPooling
    from .random_pool import RandomPooling
except ImportError:
    EdgePooling = None
    LouvainPooling = None
    RandomPooling = None
