#!/bin/bash

# Run as: 
# CS_download.sh guest your_email

# check if two files were passed as argument
if (( $# < 1 )); then
  echo "Please run as: CS_download your_email"
  exit
  
fi

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/ChIPseeqer/trunk ChIPseeqer/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/PAGE/trunk PAGE/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/FIRE/trunk FIRE/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/BIO-C/trunk BIO-C/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/MYSCANACE/trunk MYSCANACE/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/HEATMAP/trunk HEATMAP/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/PERL_MODULES/trunk PERL_MODULES/

svn checkout --username guest --password $1 https://pbtech-vc.med.cornell.edu/public/svn/elementolab/PERL_SCRIPTS/trunk PERL_SCRIPTS/