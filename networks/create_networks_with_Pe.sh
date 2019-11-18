#!/bin/bash

network_file=$1
network_topology=$2
Pe_value=$3
python3 $network_file $network_topology $Pe_value &
