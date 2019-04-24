#!/bin/bash
svn mkdir https://pbtech-vc.med.cornell.edu/public/svn/elementolab/$1 -m "created $1 dir"
svn mkdir https://pbtech-vc.med.cornell.edu/public/svn/elementolab/$1/trunk -m "created $1/trunk dir"
svn import $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/$1/trunk -m "imported $1 to trunk dir "
mv $1 $1-ORIGINAL
svn co https://pbtech-vc.med.cornell.edu/public/svn/elementolab/$1/trunk $1/



