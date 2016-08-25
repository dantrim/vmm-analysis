#!/bin bash

mkdir vector_lib/

cd include/

rootcint -f aDict.cxx -c a.h LinkDef.h

g++ -o libMylib.so aDict.cxx `root-config --cflags --libs` -shared -fPIC

mv *.so ../vector_lib/
mv *.pcm ../vector_lib/

#mkdir ../build/objects
#g++ -o ../build/objects/libMylib.so aDict.cxx `root-config --cflags --libs` -shared -fPIC
#cd ../build
#ln -s objects/libMylib.so .

#if [ "${is_mac_}" == "--mac" ]
#then
#    qmake -spec macx-g++ -o Makefile vmmdcs.pro
#    #qmake -spec macx-g++ -o Makefile vmmall.pro
#else
#    qmake -o Makefile vmmdcs.pro
#    #qmake -o Makefile vmmall.pro
#fi


#qmake -o Makefile vmmall.pro
