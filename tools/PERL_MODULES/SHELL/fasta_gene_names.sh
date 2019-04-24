#!/bin/bash
grep ">" $1 | sed "s/>//g"
