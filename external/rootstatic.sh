
#https://root.cern.ch/root/roottalk/roottalk00/2730.html

/bin/rm roota
/bin/rm libRoot.a
ar rv libRoot.a \
$ROOTSYS/base/src/*.o \
$ROOTSYS/cint/src/*.o \
$ROOTSYS/clib/src/*.o \
$ROOTSYS/cont/src/*.o \
$ROOTSYS/g3d/src/*.o \
$ROOTSYS/gpad/src/*.o \
$ROOTSYS/graf/src/*.o \
$ROOTSYS/gui/src/*.o \
$ROOTSYS/hist/src/*.o \
$ROOTSYS/histpainter/src/*.o \
$ROOTSYS/html/src/*.o \
$ROOTSYS/matrix/src/*.o \
$ROOTSYS/meta/src/*.o \
$ROOTSYS/minuit/src/*.o \
$ROOTSYS/net/src/*.o \
$ROOTSYS/physics/src/*.o \
$ROOTSYS/postscript/src/*.o \
$ROOTSYS/proof/src/*.o \
$ROOTSYS/rint/src/*.o \
$ROOTSYS/tree/src/*.o \
$ROOTSYS/treeplayer/src/*.o \
$ROOTSYS/treeviewer/src/*.o \
$ROOTSYS/unix/src/*.o \
$ROOTSYS/x11/src/*.o \
$ROOTSYS/x11ttf/src/*.o \
$ROOTSYS/x3d/src/*.o \
$ROOTSYS/zip/src/*.o
echo 'int  G__globalsetup() {}' >globalsetup.c
cc -c globalsetup.c
g++ -o roota \
$ROOTSYS/main/src/rmain.o \
globalsetup.o \
$ROOTSYS/gui/src/G*.o \
$ROOTSYS/html/src/G*.o \
$ROOTSYS/histpainter/src/G*.o \
$ROOTSYS/matrix/src/G*.o \
$ROOTSYS/treeplayer/src/G*.o \
$ROOTSYS/treeviewer/src/G*.o \
$ROOTSYS/x11/src/G*.o \
$ROOTSYS/x3d/src/G*.o \
$ROOTSYS/postscript/src/G*.o \
libRoot.a \
$HOME/xpm-3.4j/lib/libXpm.a \
/usr/X11R6/lib/libX11.a -lm -ldl -static
