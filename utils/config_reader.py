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

INT = re.compile('^[-+]?\d+$')
FLOAT = re.compile('^[-+]?\d*\.\d*$')
SYMBOL = re.compile('^\D$')
STRING = re.compile('[\w\.\/]+')
STRING2 = re.compile("\'[\w\.\/]+\'")
LIST = re.compile('^\[.+?\]$')

# value:string, recognize value type by rx match
def convert(val):
    if INT.match(val) is not None:
        return int(val)
    if FLOAT.match(val) is not None:
        return float(val)
    if SYMBOL.match(val) is not None:
        return val
    if STRING.match(val) is not None:
        return val
    if STRING2.match(val) is not None:
        return val[1:-1]  # rm apostrophes
    if LIST.match(val) is not None:
        items = [item.strip() for item in val[val.find('[')+1:val.find(']')].split(',')]
        return [convert(item) for item in items]
    else:
        print 'Warning, unrecognized type of value: "', val, '"'
        return None

class Config(object):
    def __init__(self, path_to_config):
        if os.path.exists(path_to_config) is False:
            print "Give correct path/name to configuration file"
            sys.exit(-1)
        self.var = {}
        self.order = []
        with open(path_to_config) as f:
            print '---------------- Configuration ----------------'
            lines = [line.strip() for line in f.readlines() if line.startswith('#') is False]
            for pos, line in enumerate(lines):
                vv = [v.strip() for v in line.split('=')]
                if len(vv) != 2:
                    logging.error("Unexpected row format in configuration file: " + line)
                    sys.exit(-1)
                if vv[1].find('<') > -1:
                    logging.warning("Variable not set: " + vv)
                    continue
                self.var[vv[0]] = convert(vv[1])
                self.order.append((vv[0], pos))
        logging.info(self.var)
