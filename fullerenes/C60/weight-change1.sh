#!/bin/bash
awk 'END{print $2,$3,$4,$5,$6,$7,$8,$9}'  tau1 > file1
xargs -n1 < file1 > last-weights1

