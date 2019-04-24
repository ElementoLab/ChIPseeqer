#!/bin/bash
mysql FLYDB -e "select * from FLY where ID = '$1' or CGID = '$1' or NAME = '$1'\G"
