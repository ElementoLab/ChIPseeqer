#!/bin/bash
aa=`perl ~/PERL_MODULES/SCRIPTS/revcomp.pl $1`
egrep "($1|$aa)" $2

