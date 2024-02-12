"""
Recipe to filters a field from a plotfile
"""

class Colander():
    
    def __init__(self, key):
        self.key = key

    def __call__(self, data):
        return data[self.key]
        
