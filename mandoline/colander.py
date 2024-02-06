"""
Recipe to filters a field from a plotfile
"""

def colander_generator(key):
    
    def colander(data):
        return data[key]
    colander.__doc__ = key
    
    return colander
        