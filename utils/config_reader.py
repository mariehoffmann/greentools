'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================

    Task: Read a configuration file and store settings.
    
    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================

'''

import os
import sys
import logger

class Config(object):
    self __init__(self, path_to_config):
        if os.path.exists(path_to_config) is False:
            print "Give correct path/name to configuration file"
            sys.exit(-1)
        with open(path_to_config) as f:
            lines = [line for line in f.readlines() if line.startswith('#') is False]
            for line in lines:
                vv = line.split('=')
                if len(vv) != 2:
                    print "Unexpected row format in configuration file: ", line
                    sys.exit(-1)
                # todo: how convert text into variables
                if vv[1].startswith('<'):
                    print "Warning: ", vv[0], " is not set"
                    continue
                self.vv[0] = vv[1]



