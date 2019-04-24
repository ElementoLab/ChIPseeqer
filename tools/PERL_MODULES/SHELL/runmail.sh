#!/bin/bash
echo $*
mail elemento@princeton.edu -s "$* @ $HOSTNAME" <<EOF
Machine : $HOSTNAME
Job     : $*
EOF

