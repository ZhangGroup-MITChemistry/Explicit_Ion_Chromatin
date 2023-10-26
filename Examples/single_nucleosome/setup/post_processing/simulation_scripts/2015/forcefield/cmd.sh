#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Sep 15 17:55:32 2022
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash

# $i 35 lj/cut 0.239 sigma_ij=(sigma_ions+5.7)/2 sigma_ij*(2)^(1/6) 
for ((i=14; i<=34; i++)); do echo "$i 35 lj/cut 0.239 4.097 4.5987" >> tmp.txt; done
for ((i=14; i<=34; i++)); do echo "$i 36 lj/cut 0.239 4.85 5.4439" >> tmp.txt; done
for ((i=14; i<=34; i++)); do echo "$i 37 lj/cut 0.239 5.089 5.7122" >> tmp.txt; done
for ((i=14; i<=34; i++)); do echo "$i 38 lj/cut 0.239 4.5 5.0511" >> tmp.txt; done
