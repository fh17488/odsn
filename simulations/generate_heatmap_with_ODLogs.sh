#!/bin/bash

java -classpath . RAwNetworkModel "ring-lattice" 0.012 1.96
java -classpath . RAwNetworkModel "small-world" 0.012 1.96
java -classpath . RAwNetworkModel "random" 0.012 1.96

