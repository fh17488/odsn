#!/bin/bash

network_topology=$1
key_input_params=$2
counter=1
for Pe in '0' '0.006' '0.012' '0.018' '0.024' '0.03' '0.036' '0.042' '0.048' '0.054' '0.06' '0.066' '0.072' '0.078' '0.084' '0.09' '0.096' '0.102' '0.108' '0.114' '0.12' '0.126' '0.132' '0.138' '0.144' '0.15' '0.156' '0.162' '0.168' '0.174' '0.18' '0.186' '0.192' '0.198' '0.204' '0.21' '0.216' '0.222' '0.228' '0.234' '0.24' '0.246' '0.252' '0.258' '0.264' '0.27' '0.276' '0.282' '0.288' '0.294'
do
    for U in '0.000' '0.040' '0.080' '0.120' '0.160' '0.200' '0.240' '0.280' '0.320' '0.360' '0.400' '0.440' '0.480' '0.520' '0.560' '0.600' '0.640' '0.680' '0.720' '0.760' '0.800' '0.840' '0.880' '0.920' '0.960' '1.000' '1.040' '1.080' '1.120' '1.160' '1.200' '1.240' '1.280' '1.320' '1.360' '1.400' '1.440' '1.480' '1.520' '1.560' '1.600' '1.640' '1.680' '1.720' '1.760' '1.800' '1.840' '1.880' '1.920' '1.960' 
    do    
        java -classpath . RAwNetworkModel $network_topology $Pe $U $key_input_params &
        counter=$((counter+1))
        if (( $counter % 17 == 0 ))
        then
            sleep 6s
        fi
    done
done

