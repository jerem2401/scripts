#! /bin/bash

name="$(grep "^model name" /proc/cpuinfo | head -n1 | cut -d: -f2)"
# echo "$(basename $0): found model name = $name" >&2

if echo "$name" | egrep -q "Xeon.*E5-2650.*v4|CPU E5-2630.*v4"; then
    # Broadwell
    A=broadwell
elif echo "$name" | egrep -q 'Xeon.*E5-2670|Xeon.*E3-1270|Xeon.*E5-2650|Xeon.*4214'; then
   #  Sandy Bridge
   A=sandy-bridge
elif echo "$name" | egrep -q 'Xeon.*Gold'; then
   A=cascadelake
elif echo "$name" | egrep -q 'AMD .*Processor 6378|AMD.* 6220'; then
   # Interlagos
   A=interlagos
elif echo "$name" | egrep -q 'Xeon.* E5540'; then
   # Nehalem
   A=nehalem
elif echo "$name" | egrep -q 'AMD.*6174|AMD.*6136'; then
   # Magny-Cours
   A=magny-cours
elif echo "$name" | egrep -q "Xeon.*X5650"; then
    # GPU nodes GeForce GTX 480
    A=gtx480
elif echo "$name" | egrep -q "Xeon.*E5-2680"; then
    # Haswell 
    #  --- FOR NOW, USE THE SANDY-BRIDGE INSTALLATION
    #
    A=sandy-bridge
elif echo "$name" | egrep -q "Xeon.*E5-2650.*v4"; then
    # Broadwell
    A=broadwell
elif echo "$name" | egrep -q "AMD Opteron.*Processor 6276"; then
    # GWDG login machine
    A=gwdg-login
else
   echo "$(basename $0): Could not detect Architechture: $name" >&2; exit 1
fi

echo $A
