import utils
from typing import List
from collections import defaultdict

# ==============================================================================================================
# Node Class
# NOTE: Node class is an adaption of the TreeNode class from scikit-bio.
# https://github.com/biocore/scikit-bio/blob/0.2.3/skbio/tree/_tree.py#L53
# ==============================================================================================================
class Node:
    '''
    Stores an instance of a tree data structure.

    Parameter(s):
        name (str, default=None): name of the species or type.
        distance (float, default=None): evolutionary distance between species/nodes.
        parent (Node, default=None): parent node.
        children (List[Node], default=None): nodes of existing children.

    Output(s): None
    '''
    def __init__(self, name:str=None, distance:float=None, parent:'Node'=None, children:List['Node']=None):
        self.name = name
        self.distance = distance
        self.parent = parent

        self._tip_cache = {}
        self._non_tip_cache = {}
        self._registered_caches = set()
        self.children = []

        if children is not None:
            self.extend(children)

    # ----------------------------------------------------------------------------------------------------------
    def __str__(self):
        '''
        Adds children to the tree.
        
        Parameter(s): None
            
        Output(s):
            A string of the tree data comprised of names and distances
        '''
        return self.name
    
    # ----------------------------------------------------------------------------------------------------------
    def is_root(self):
        '''
        Returns True if the current node is the root.
        
        Parameter(s): None
        
        Output(s):
            True if the current node is the root, else False if the node is not the root.
        '''

        return self.parent is None
    
    # ----------------------------------------------------------------------------------------------------------
    def root(self):
        '''
        Gets the root of the tree.
        
        Parameter(s): None
        
        Output(s): None
        '''

        node = self
        while not node.is_root():
            node = node.parent

        return node
    
    # ----------------------------------------------------------------------------------------------------------
    def invalidate_caches(self, attr:bool=True):
        
        if not self.is_root():
            self.root().invalidate_caches()
        else:
            self._tip_cache = {}
            self._non_tip_cache = {}

            if self._registered_caches and attr:
                for n in self.traverse():
                    for cache in self._registered_caches:
                        if hasattr(n, cache):
                            delattr(n, cache)

    # ----------------------------------------------------------------------------------------------------------
    def create_caches(self):
        '''
        Constructs a lookup cache for node names.
        
        Parameter(s): None
        
        Output(s): None
        '''

        if not self.is_root():
            self.root().create_caches()
        else:
            if self._tip_cache and self._non_tip_cache:
                return

            self.invalidate_caches(attr=False)

            tip_cache = {}
            non_tip_cache = defaultdict(list)

            for node in self.postorder():
                name = node.name

                if name is None:
                    continue

                if node.is_tip():
                    if name in tip_cache:
                        raise utils.InvalidInput(f"Tip with name {name} already exists!")

                    tip_cache[name] = node
                else:
                    non_tip_cache[name].append(node)

            self._tip_cache = tip_cache
            self._non_tip_cache = non_tip_cache

    # ----------------------------------------------------------------------------------------------------------
    def _adopt(self, node:'Node'):
        self.invalidate_caches()

        if node.parent is not None:
            node.parent.remove(node)
        node.parent = self

        return node
    # ----------------------------------------------------------------------------------------------------------
    def extend(self, nodes:List['Node']):
        '''
        Adds children to the tree.
        
        Parameter(s):
            nodes (List[Node]): a list of children to be added to the tree.
            
        Output(s): None
        '''

        self.children.extend([self._adopt(n) for n in nodes])
    
    # ----------------------------------------------------------------------------------------------------------
    def find(self, name:str):
        '''
        Searches for the first instance of a node with the same name.
        
        Parameter(s):
            name (str): name of the node being searched.
            
        Output(s): None
        '''
        root = self.root()

        # if what is being passed in looks like a node, just return it
        if isinstance(name, root.__class__):
            return name

        root.create_caches()
        node = root._tip_cache.get(name, None)

        if node is None:
            node = root._non_tip_cache.get(name, [None])[0]

        if node is None:
            raise utils.InvalidInput(f"Node {name} is not in self!")
        else:
            return node
    
    # ----------------------------------------------------------------------------------------------------------
    def append(self, node:'Node'):
        '''
        Appends a node to children.

        Parameter(s):
            node (Node): node being appended to list of children
        
        Output(s): None
        '''
        self.children.append(self._adopt(node))

    # ----------------------------------------------------------------------------------------------------------
    def pop(self, index:int=-1):
        '''
        Remove a Node from the tree.

        Parameter(s):
            idx (int): index location of the node being removed.
        
        Output(s):
            The node removed from the tree. 
        '''
        return self._remove_node(index)

    # ----------------------------------------------------------------------------------------------------------
    def _remove_node(self, index:int):
        '''
        Performs node removal.
        
        Parameter(s):
            idx (int): index location of the node being removed.
        
        Output(s):
            The node removed from the tree.        
        '''
        self.invalidate_caches()
        node = self.children.pop(index)
        node.parent = None

        return node

    # ----------------------------------------------------------------------------------------------------------
    def remove(self, node:'Node'):
        '''
        Removes a node from the tree.

        Parameter(s):
            node (Node): the current node being removed.

        Output(s):
            True if the node was removes, else returns false.
        '''
        for (i, curr_node) in enumerate(self.children):
            if curr_node is node:
                self._remove_node(i)
                return True
            
        return False
    
    # ----------------------------------------------------------------------------------------------------------
    def traverse(self, self_before:bool=True, self_after:bool=False, include_self:bool=True):
        '''
        Iterates over descendants and returns them.
        
        Paramter(s):
            self_before (bool, default=True): includes each node before its descendants if True.
            self_after (bool, default=False): includes each node after its descendants if True.
            include_self (bool, default=False): include the initial node if True.
            
        Output(s):
            Yields successive Node objects.
        '''
        if self_before:
            if self_after:
                return self.pre_and_postorder(include_self=include_self)
            else:
                return self.preorder(include_self=include_self)
        else:
            if self_after:
                return self.postorder(include_self=include_self)
            else:
                return self.tips(include_self=include_self)

    # ----------------------------------------------------------------------------------------------------------
    def preorder(self, include_self:bool=True):
        '''
        Performs preorder iteration over the tree.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.
            
        Output(s):
            Yields successive Node objects.
        '''
        
        stack = [self]
        while stack:
            curr = stack.pop()
            if include_self or (curr is not self):
                yield curr
            if curr.children:
                stack.extend(curr.children[::-1])
    
    # ----------------------------------------------------------------------------------------------------------
    def postorder(self, include_self:bool=True):
        '''
        Performs postorder iteration over tree.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.

        Output(s):
            Yields successive Node objects.
        '''
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        curr_children_len = len(curr_children)
        while 1:
            curr_index = child_index_stack[-1]
            # if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                curr_children_len = len(curr_children)
                child_index_stack.pop()
                child_index_stack[-1] += 1

    # ----------------------------------------------------------------------------------------------------------
    def pre_and_postorder(self, include_self:bool=True):
        '''
        Performs iteration over tree, visiting node before and after.
        
        Parameter(s):
            include_self (bool, default=True): include the initial node if True.

        Output(s):
            Yields successive Node objects.
        '''
        # handle simple case first
        if not self.children:
            if include_self:
                yield self
            raise StopIteration
        child_index_stack = [0]
        curr = self
        curr_children = self.children
        while 1:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            # if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                # if the current child has children, go there
                if curr_child.children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.children
                    curr_index = 0
                # otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            # if there are no children left, return self, and move to
            # self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.parent
                curr_children = curr.children
                child_index_stack.pop()
                child_index_stack[-1] += 1

    # ----------------------------------------------------------------------------------------------------------
    def is_tip(self):
        '''
        Finds if the node is a tip (has no children) or not.

        Parameter(s): None

        Output(s):
            Returns True if the node is a tip and has no children, else False.
        '''
        return not self.children
    
    # ----------------------------------------------------------------------------------------------------------
    def tips(self, include_self:bool=False):
        '''
        Iterates over tips descended from self.

        Parameter(s):
            include_self (bool, default=False): include the initial node if True

        Output(s):
            Yeilds successive Node objects.
        '''
        for node in self.postorder(include_self):
            if node.is_tip():
                yield node

    # ----------------------------------------------------------------------------------------------------------
    def non_tips(self, include_self:bool=False):
        '''
        Iterates over nontips descended from self.

        Parameter(s):
            include_self (bool, default=False): include the initial node if True

        Output(s):
            Yeilds successive Node objects.
        '''
        for node in self.postorder(include_self):
            if not node.is_tip():
                yield node