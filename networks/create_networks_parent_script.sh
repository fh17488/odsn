#!/bin/bash

sleep 60m
./create_networks.sh "ring-lattice" ./create_network_topology_wExtremists_K16_R0.py
sleep 60m
./create_networks.sh "ring-lattice" ./create_network_topology_wExtremists_K16_R0_25.py
sleep 60m
./create_networks.sh "ring-lattice" ./create_network_topology_wExtremists_K16_R0_75.py
sleep 60m
./create_networks.sh "ring-lattice" ./create_network_topology_wExtremists_K16_R1.py
sleep 60m
./create_networks.sh "small-world" ./create_network_topology_wExtremists_K128.py

