#!/usr/bin/python3

import threading
import time
import numpy as np
import logging
import random as rg
import sys
from timeit import default_timer as timer
import csv
import os


# argv[1]: Network Topology
# argv[2]: Specific Data Directory
# argv[3]: Key Input Parameters
def main():
    topology = topology = str(sys.argv[1])
    logging.basicConfig(level=logging.INFO)
    title = topology+"_heatmap_N1000_"+str(sys.argv[3])
    data_directory = "../Results/"+str(sys.argv[2])+"/"
    outfile = title+".csv"
    outfile_handle = open(outfile, 'w')
    outstr = ""
    u = ['0.0','0.04','0.08','0.12','0.16','0.2','0.24','0.28','0.32','0.36','0.4','0.44','0.48','0.52','0.56','0.6','0.64','0.68','0.72','0.76','0.8','0.84','0.88','0.92','0.96','1.0','1.04','1.08','1.12','1.16','1.2','1.24','1.28','1.32','1.36','1.4','1.44','1.48','1.52','1.56','1.6','1.64','1.68','1.72','1.76','1.8','1.84','1.88','1.92','1.96']
    pe = ['0.0','0.006','0.012','0.018','0.024','0.03','0.036','0.042','0.048','0.054','0.06','0.066','0.072','0.078','0.084','0.09','0.096','0.102','0.108','0.114','0.12','0.126','0.132','0.138','0.144','0.15','0.156','0.162','0.168','0.174','0.18','0.186','0.192','0.198','0.204','0.21','0.216','0.222','0.228','0.234','0.24','0.246','0.252','0.258','0.264','0.27','0.276','0.282','0.288','0.294']
    u_steps = 50
    pe_steps = 50
    for pe_i in range(pe_steps):
        for u_i in range(u_steps):
            pe_val = pe[pe_i]
            u_val = u[u_i]
            infile = data_directory + topology + "_network_" + str(pe_val) + "_"  + str(u_val) +  "_MetricYAvgLog" + str(sys.argv[3]) +".txt"
            print(infile)
            infile_handle = open(infile, 'r')
            if u_i < (u_steps-1):
                outstr += infile_handle.read().rstrip('\n') + ", "
            else:
                outstr += infile_handle.read().rstrip('\n') + "\n"
            infile_handle.close()
    outfile_handle.write(outstr)
    outfile_handle.close()
    return


main()
sys.exit(0)

