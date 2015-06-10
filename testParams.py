__author__ = 'pbmanis'


class Params(object):
    """
    utility class to create parameter lists...
    create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
    reference like p.abc, p.defg, etc.
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def getkeys(self):
        """
        Get the keys in the current dictionary
        """
        return(self.__dict__.keys())

    def haskey(self, key):
        """
        Find out if the param list has a specific key in it
        """
        if key in self.__dict__.keys():
            return True
        else:
            return False
    def todict(self):
        """
        convert param list to standard dictionary
        Useful when writing the data
        """
        r = {}
        for dictelement in self.__dict__:
            if isinstance(self.__dict__[dictelement], Params):
                #print 'nested: ', dictelement
                r[dictelement] = self.__dict__[dictelement].todict()
            else:
                r[dictelement] = self.__dict__[dictelement]
        return r

    def show(self):
        """
        print the parameter block created in Parameter Init
        """
        print "--------    Parameter Block    ----------"
        for key in self.__dict__.keys():
            print "%15s = " % (key), eval('self.%s' % key)
        print "-------- ---------------------- ----------"

z= Params(m=0,n=1,p=2)
y = Params(x=9,y=3,z=z)
x = Params(a=1, y=y)
print 'x: ', x
u=x.todict()
print 'u: ', u
x.show()