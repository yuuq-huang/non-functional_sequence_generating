class TrieNode:
 
    def __init__(self, char):
        self.char = char
        self.is_end = False
        self.children = {}
 
class Trie(object):
 
    def __init__(self):
        self.root = TrieNode("")
     
    def insert(self, word):
        node = self.root
        
        for char in word:
            if char in node.children:
                node = node.children[char]
            else:
                new_node = TrieNode(char)
                node.children[char] = new_node
                node = new_node
         
        node.is_end = True
         
    def check_if_end(self, x):
        
        node = self.root
        
        if not node.children:
            return False
         
            
        for char in x:
            if char in node.children:
                node = node.children[char]
            else:
                return False
       
        if node.is_end:
           return True
        else:
           return False
    
    def get_node(self,x):
        node = self.root
        
        for char in x:
            if char in node.children:
                 node = node.children[char]
            else:
                return []
            
        self.output=[]
        self.output.append(node.char)
        return self.output
    
    def search(self,x):
        node = self.root
        
        for char in x:
            if char in node.children:
                node = node.children[char]
            else:
                return False
        if self.check_if_end(x):
            return True
        else:
            return False
    def search_children(self,x):
        node = self.root
        
        for char in x:
            if char in node.children:
                node = node.children[char]
            else:
                return []
        #self.output=[]
        #self.output.append(node.children)
        return node.children
