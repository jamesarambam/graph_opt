"""
++++++++++++++++++++++++++++++++++
Author : James Arambam
Date   : 28 Jan 2016

++++++++++++++++++++++++++++++++++
"""
import os
import auxLib as ax

# =============================================================================== #

f = "49_4_2_1.wcsp"
os.system("python graphOpt.py data/"+f+" -min ")
exit()

file = ax.listFiles("data/")
for f in file:
    os.system("python graphOpt.py data/"+f+" -min ")
    exit()
exit()
