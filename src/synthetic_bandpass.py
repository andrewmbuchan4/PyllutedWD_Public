#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

class Bandpass(Enum):
    u = 0
    g = 1
    r = 2
    i = 3
    z = 4
    
    def __str__(self):
        return self.name
