import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from time import time,ctime


plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 15}

plt.rc('font', **font)

print ' Python plotting...'

pos=[]
ins=[]
nins=[]

file_name = sys.argv[1]
chromosome = sys.argv[2]
fp = open(file_name, 'r')
for line in fp:
    col= line.split()
    pos.append(float(col[1]))
    ins.append(float(col[2]))
    nins.append(float(col[3]))
fp.close()

x = np.array(pos)
y = np.array(ins)
ny= np.array(nins)

myfig = PdfPages(chromosome + '.pdf') 

fig = plt.figure(figsize=(13,6.5))

# change hor scale to include all genome length?

fig.add_subplot(1,2,1)
plt.ylabel('Insert Size')
plt.xlabel('Position on Chromosome '+ chromosome)
#text(x, y, s, fontsize=12)
plt.plot(x,y,'ro')
plt.axhline(y=0) 
plt.tight_layout()
#plt.show()

fig.add_subplot(1,2,2)
plt.ylim((-1,1))
plt.ylabel('Normalized Insert Size')
plt.xlabel('Position on Chromosome '+ chromosome)
#text(x, y, s, fontsize=12)
plt.plot(x,ny,'ro')
plt.axhline(y=0)
plt.tight_layout()

plt.show()


#myfig = PdfPages(chromosome + '.pdf') 
myfig.savefig(fig)
myfig.close()

