exec 2>&1
#!/bin/bash
# Usage:     adu_sequence row column
# Dumps the corresponding adu sequence 

h5dump -d "/10ms/stack[$1,$2,0;;1,1,100;]" -A 0 CCD_calibration.hdf5 | \
    grep "([0-9][0-9],[0-9][0-9],\(0\|[0-9][0-9]\)):" | \
    sed -e "s/[[:space:]]*([0-9][0-9],[0-9][0-9],\(0\|[0-9][0-9]\)):[[:space:]]*//g" |\
    sed -e "s/,$//g" | sed -e "s/,[[:space:]]*/\n/g"
:
