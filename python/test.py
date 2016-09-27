#!/usr/bin/env python

class Test :
    def __init__(self, no) :
        self.number = no

    def Print(self) :
        print self.number

    def __eq__(self, rhs) :
        return self.number == rhs.number



def main() :
    print "TESTING"

    tmp = []
    x = Test(0)
    y = Test(2)
    z = Test(1)
    q = Test(43)

    tmp.append(y)
    tmp.append(x)
    tmp.append(q)
    tmp.append(z)

    b = Test(0)
    if b==x :
        print "b equals x"
    print "b number %d"%b.number
    print "x number %d"%x.number

    tmp.sort(key=lambda x: x.number, reverse=False)

    for hit in tmp :
        hit.Print()

    blah = Test(2)
    blah2 = Test(3)
    yep = (abs(blah.number-blah2.number)==1)
    print "yep: ", yep

if __name__ == "__main__" :
    main()
