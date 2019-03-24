import os

f = open('/home/rabrahm/hd72673.list','r')
lines = f.readlines()

for line in lines:
    night = line.split()[0]
    os.system('python ferospipe.py ' + night + ' -npools 10 -do_class -o2do HD72673 -avoid_plot')
