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
import re
import sys
import logging

NUMBER = re.compile('[-+]?\d*\.?\d+')
STRING = re.compile('[\w\.\/]+')
LIST = re.compile('\[.+?\]')

# value:string, recognize value type by rx match
def convert(val):
    print 'convert called with val = ' + val
    if NUMBER.match(val) is not None:
        return float(val)
    if STRING.match(val) is not None:
        return val
    if LIST.match(val) is not None:
        items = [item.strip() for item in val[val.find('[')+1:val.find(']')].split(',')]
        print items
        return [convert(item) for item in items]

class Config(object):
    def __init__(self, path_to_config):
        if os.path.exists(path_to_config) is False:
            print "Give correct path/name to configuration file"
            sys.exit(-1)
        self.var = {}
        with open(path_to_config) as f:
            lines = [line for line in f.readlines() if line.startswith('#') is False]
            for line in lines:
                vv = line.split('=')
                if len(vv) != 2:
                    logging.error("Unexpected row format in configuration file: " + line)
                    sys.exit(-1)
                # todo: how convert text into variables
                if vv[1].find('<') > -1:
                    logging.warning("Variable not set: " + vv)
                    continue
                self.var[vv[0]] = convert(vv[1])
        logging.info(self.var)
