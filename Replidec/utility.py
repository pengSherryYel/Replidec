#!/usr/bin/env python
# coding: utf-8
# authors: sherry peng, torben sanders, erfan khamespanah
# mail: xue.peng@helmholtz-muenchen.de
# date: 2026.07.09

import os
from subprocess import Popen, PIPE
import sys


def mkdirs(dirname):
    '''Safely creates output directories if they do not already exist.'''
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def checkEnv(sft):
    '''
    Validates that required third-party binaries are accessible in the system PATH.

    Uses standard system commands to check for binary availability. If a dependency
    cannot be resolved, it raises a clean critical alert and exits to prevent
    pipeline failures down the line.
    '''
    cmd = "which %s" % sft

    # Direct stdout and stderr streams to execution PIPEs to preserve clean console logging
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    obj.wait()

    # An execution return code other than 0 indicates missing dependencies
    if obj.returncode != 0:
        print(f"\n[CRITICAL] Dependency Missing: '{sft}' could not be resolved in your current environment.")
        print(f"Please ensure '{sft}' is installed and available in your PATH variable.\n")
        sys.exit(1)