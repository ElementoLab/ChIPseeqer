#!/bin/sh
fix_random_sequence_data.pl $1 > /tmp/tmpfix
mv /tmp/tmpfix $1
