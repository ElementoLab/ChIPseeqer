#!/bin/bash

d=`cygpath.exe -d $CHIPSEEQERDIR/TREEVIEW/TreeView.jar`
p=`cygpath.exe -d $CHIPSEEQERDIR/$1`

java -jar -Xmx800m $d -r $p
