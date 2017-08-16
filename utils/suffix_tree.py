'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================

    Compressed suffix tree implementation with helper structures for e.g. suffix
    links of alltogether size |CSA| + 6n + o(n). These structures are needed to
    answer questions about longest common substrings of two strings, matching
    statistics, or to enumerate maximal repeated substrings efficiently.
    The compressed suffix tree and its construction are described here:
    K. Sadakane: "Compressed suffix trees with full functionality", 2007.

    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================
'''

# used by SuffixTree
class SuffixArray(object):
    def __init__(self):
        pass

class SuffixTree(object):

    # inner node class
    class Node(object):
        def __init__(self, parent=None, label=[-1,-1], children=[]):
            self.parent = parent      # ptr to outer node array
            self.children = children  # list of pointers to node array
            # label on edge from ancestor to this node stored as text ptrs
            self.label = label
            self.sl = None            # suffix link to another node


        def label_len(self):
            return label[1] - label[0]

    def __init__(self, text=''):
        self.__text = text  # encoded in nlog|A| bits
        self.__root = 0  # starting point for tree topology
        self.__nodes = [self.Node()] # stored in inorder ordering
        self.__node_depths = []  # store string depth for each internal node
        self.__edge_labels = []
        self.__suffix_links = []

    # private helper functions to resolve indices to node objects
    def __children(self, v_idx):
        return [self.__nodes[w_idx] for w_idx in self.nodes[v_idx].children]

    def __label(self, v_idx):
        return self.__text[self.nodes[v_idx].label[0]:self.nodes[v_idx].label[0]]

    def __parent(self, v_idx):
        if self.__root == v_idx:
            return None
        return self.__nodes[self.__nodes[v_idx].parent]

    def __root(self):
        return self.__nodes[self.__root]

    def __sl(self, v_idx):
        return self.__nodes[v_idx].sl

    # public member functions are dealing internally with node indices
    # meta information is currently stored in the node objects (in __nodes)
    def root(self):
        return self.__root()

    # output True if node v is a leaf
    def isleaf(self, v):
        if len(self.__children(v)) == 0:
            return True
        return False

    # return child as node_list index of v if its edge label starts with character c
    def child(self, v, c):
        for u in self.__children(v):
            if self.__label(u).startswith('c'):
                return self.__nodes[u]
        return None

    # return first child of node v
    def firstchild(self, v):
        if len(self.__children(v)) > 0:
            return self.__children(v)[0]
        return None

    # return next sibling of node v
    def sibling(self, v):
        if self.__root == v:
            return None
        for i, u in enumerate(__children(__parent(v))[:-1]):
            if u == v:
                return self.__nodes[self.__parent(v).children[i+1]]
        return None

    # return parent of node v
    def parent(self, v):
        return self.__parent(v)

    # return dth character of edge label pointing to v
    def edge(v, d):
        if __debug__:
            if v == self.root: raise AssertionError("Root node has no label.")
            if d >= self.__label(v).label_len(): raise AssertionError("Index d exceeds node label length.")
        return self.__label(v)[d]

    # return string depth (total string length) of node v
    def depth(v):
        d = len(self.__label(v))
        current = v.parent
        while current != self.root:
            d += len(self.__label(current))
            current = current.parent
        return d

    # return lowest common ancestor between nodes v and w
    def lca(v, w):
        anc_v = set([v])
        current = v
        while current != self.root:
            current = self.__nodes[current].parent
            anc_v.add(current)
        current = w
        while current != self.root:
            current = self.__nodes[current].parent
            if current in anc_v:
                return current
        return self.root

    # return node w that is pointed to by suffix link from node v
    def sl(v):
        return self.__sl(v)

    # compute node depths of final tree for each internal node with bfs
    def compute_node_depths(self):
        if __debug__:
            if self.root == None: raise AssertionError("Tree is empty.")
        self.node_depths = []
        stack = [(v, v.label_len()) for v in self.__children(self.root)]
        while len(stack) > 0:
            current = stack.pop()
            if not self.isleaf(current[0]):
                for u in self.__nodes[current[0]].children:
                    self.__node_depths.append((u, current[1] + self.__label(u)))

    # build index ms of length 2|t|, append for each i the binary string
    def compressed_suffix_tree(t):
        pass

# tests
if __name__ == '__main__':
    st = SuffixTree()
