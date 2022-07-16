#!/usr/bin/env python
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6

import os
from subprocess import Popen


def mkdirs(dirname):
    '''
    makedirs
    '''
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def checkEnv(sft):
    '''
    check software in the PATH
    '''
    cmd = "which %s" % sft
    status = Popen(cmd, shell=True)
    status.wait()
    if status:
        print("%s exist" % sft)
    else:
        print("Please add %s in your PATH" % sft)
