import os
import time
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Input utilisateur
# -----------------
lh = 0.5
lv = 0.5
H = 6.0
E = 210e9
nu = 0.2
t = 0.05
# lst_coins contient des 2-tuples.
lst_coins = [(13.5/2,6),(10,5.5/2),(10,-5.5/2),(13.5/2,-6),(-13.5/2,-6),(-10,-5.5/2),(-10,5.5/2),(-13.5/2,6)]
recalculate = True

# Traitement
# ----------

fileIn = open('input_batardeau.txt','w')

fileIn.write(str(lh)+"\n")
fileIn.write(str(lv)+"\n")
fileIn.write(str(H)+"\n")
fileIn.write(str(E)+"\n")
fileIn.write(str(nu)+"\n")
fileIn.write(str(t)+"\n")
fileIn.write(str(len(lst_coins))+"\n")
for k in range(len(lst_coins)):
    fileIn.write(str(lst_coins[k][0])+"\n")
    fileIn.write(str(lst_coins[k][1])+"\n")

fileIn.close()
#'''
if recalculate:
    try:
        os.remove('output_batardeau.csv')
    except:
        "Rien"
    t0 = time.time()
    a = os.system(os.getcwd()+r'/batardeau_engine')
    print('[s] '+str(round(time.time()-t0,3)))


df = pd.read_csv('output_batardeau.csv')

max_coord = np.max(np.abs(df[['x','y','z']]))
max_displ = np.max(np.abs(df[['ux','uy','uz']]))

fig4 = plt.figure()
ax4 = fig4.add_subplot(111, projection='3d')
fact = 0.2*max_coord/max_displ
scat4 = ax4.scatter([x+fact*ux for x,ux in zip(df['x'].tolist(),df['ux'].tolist())],
                    [y+fact*uy for y,uy in zip(df['y'].tolist(),df['uy'].tolist())],
                    [z+fact*uz for z,uz in zip(df['z'].tolist(),df['uz'].tolist())],
                    label=r'Déformée ($\times$'+str(int(fact))+')')
plt.axis('equal')
plt.legend()
ax4.set_xlabel('X [m]')
ax4.set_ylabel('Y [m]')
ax4.set_zlabel('Z [m]')
plt.show()
#'''

