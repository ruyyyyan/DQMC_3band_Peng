import commands
import shutil
import os
import sys
import time
import linecache
import math

for idx in range(26455,26459+1):
    cmd = "qdel " + str(idx)
    os.system(cmd) 
