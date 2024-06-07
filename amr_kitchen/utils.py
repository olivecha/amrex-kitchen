class TastesBadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# TODO: replace with more specific errors
# depending on where the problem is
class BadTastingHeadersError(Exception):
    pass

class BadTastingBinariesError(Exception):
    pass                    

