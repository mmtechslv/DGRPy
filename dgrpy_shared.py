#All of the function/method parameters must begin with i(stands for input) and must follow with capital case. For example, myfunction(iParemeter1,iParameter2,etc)
#Following is the stelodemOTU class that perfom loading of OTU file, manipulating it, data filtering etc.


class dgrpy():
    def __init__(self,iFilename,iSeparator=',',iTaxaColumn=False):
        self.filename = iFilename