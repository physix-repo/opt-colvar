#!/bin/bash
awk 'END{print $2,$3,$4,$5,$6,$7,$8,$9}'  tau > file
xargs -n1 < file > last-weights

