'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================

    Set of helper functions to predict chemical primer pair properties.

    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================
'''

import math

# Wallace rule to compute primer melting temperature
def primer_melt_wallace(primer):
    return 2*(sum([1 for nucl in primer if nucl in ['a', 't']])) + 4*(sum([1 for nucl in primer if nucl in ['c', 'g']]))

# salt-adjusted method to compute primer melting temperature
# input primer:string sequence, Na:float molar Natrium ion concentration
def primer_melt_salt(primer, Na):
    primer = primer.lower()
    return 100.5 + 41.0*(sum([1 for nucl in primer if nucl in ['c', 'g']]))/len(primer) - 820.0/len(primer) + 16.6*math.log10(Na))
