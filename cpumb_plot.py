#!/home/usersJH/jeremy/Downloads/yes/bin/python

import os
import re
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import glob
from readplumed import *
import sys
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_files with spaces inbetween, should not be used with -umb', nargs='+', action='store', dest='f')
    parser.add_argument('-umb', help='specify that umb plot mode should be turn on', action='store_true', dest='umb')
    parser.add_argument('--a', help='provide alpha if not in directory_name under certain format', action='store', dest='a', type=float)
    parser.add_argument('--o', help='give path for plot output', action='store', dest='o')
    args = parser.parse_args()


    if args.f != None:

        files=args.f

        cwd = os.getcwd()
        cwds = cwd.split('/')
        
        if args.a != None:
            alpha=float(args.a)
        else:
            alpha = re.search(r"(?<=_a).*(?=_deq)", cwds[-1]).group()
            print("careful alpha guessed from working directory name:"+str(alpha))
	        

        #function
        a=float(alpha)
        deq=0.05
        
        x = np.arange(0.001, 0.180, 0.005)
        y = np.arange(0.001, 0.180, 0.005)
        X, Y = np.meshgrid(x, y)
        
        def nCV(X,Y):
            Z = ((deq/X)**a)-((deq/Y)**a)
            return Z
        
        Z=nCV(X,Y)

        ###################################Special values################################################
        Xeq=0.001376
        Yeq=0.098296
        
        Xax=0.098361
        Yax=0.001384
        
        eqCVval=nCV(Xeq,Yeq)
        axCVval=nCV(Xax,Yax)
        #################################################################################################

        RMSDeq = []
        RMSDax = []
        phi = []
        psi = []
            
        fig=plt.figure(1, figsize=(10,6*len(files)))
        fig=plt.figure(1)

        for pltn, i in zip(range(1, len(files)*2+1, 2), files):

            RMSDeq = rdata(i)[0]
            idd=[elem for elem in range(0, len(RMSDeq),1)]
            RMSDeq = [RMSDeq[i] for i in idd]
            RMSDax = rdata(i)[1]
            RMSDax = [RMSDax[i] for i in idd]
            phi = rdata(i)[2]
            phi = [phi[i] for i in idd]
            psi = rdata(i)[3]
            psi = [psi[i] for i in idd]
            
            plt.subplot(len(files), 2, pltn)
            levels=[-2,-1,-0.5,-0.25,-0.1,-0.05,0,0.05,0.1,0.25,0.5,1,2]
            levels+=[axCVval,eqCVval]
            levels=sorted(levels)
            CS = plt.contour(X, Y, Z, levels=levels)
            plt.clabel(CS, inline=1, fontsize=10)
            plt.xlabel('RMSDeq')
            plt.ylabel('RMSDax')
            plt.scatter(RMSDeq, RMSDax, c=list(range(1, len(RMSDax)+1, 1)), cmap=plt.cm.get_cmap('rainbow'), s=10)
            plt.title('RMSDs of trajectory relative to state ax and eq \n & projection of the nCV')
            plt.xlim(0.001, 0.180)
            plt.ylim(0.001, 0.180)
            plt.gca().set_aspect('equal', adjustable='box')

            plt.subplot(len(files), 2, pltn+1)
            plt.xlim(-3.2, 3.2)
            plt.ylim(-3.2, 3.2)
            plt.title('phi & psi dihedral angles')
            plt.xlabel('phi')
            plt.ylabel('psi')
            plt.scatter(phi, psi, c=list(range(1, len(RMSDax)+1, 1)), cmap=plt.cm.get_cmap('rainbow'), s=10)
            plt.gca().set_aspect('equal', adjustable='box')

        plt.tight_layout(pad=3)
        fig.suptitle(cwds[-1])
        if args.o == None:
            plt.savefig(cwds[-1]+".jpeg")
        else :
            plt.savefig(args.o)
        plt.close()

    else:
        if args.umb == True:

            print("umb plot mode turned on")
            print(
            '''
            You've specified -umb: umb plot mode turned on

            Careful: this script generate n pdf files, with n=nbr of umb windows/15 (+1 if left over)
            Then if you have a number of windows >> 75, you should think of rearranging this script"
            '''
            )
            #path of umb working dir with all subdir called "E_*" for each um windows with diff ksi values
            #(out,err) = subprocess.Popen(['pwd'],stdout=subprocess.PIPE).communicate()
            #path=out.decode().strip()+'/' 
            #path=out.decode().strip().replace('analysis','')+"setup/"
            path=str(sys.argv[1])

            #list of all colvar file outputs from plumed, format: "colvar.*"
            #colvar_list=[os.path.join(os.path.join(path, elem), x)
            #            for elem in os.listdir(path)
            #            if os.path.isdir(os.path.join(path, elem)) == True
            #            if "E_" in elem
            #            for x in os.listdir(os.path.join(path, elem))
            #            if "colvar" in x]
            colvar_list=[os.path.join(f, x)
                    for f in glob.glob(path + "**/", recursive=True)
                    for x in os.listdir(f)
                    if re.search(r"colvar_([0-9]|-|\.)", x) !=None]

            #list of ksi values from which we've restrained the coordinates, take this values from umbsetup.sh ?
            #CVval=[elem.split('_')[1] for elem in os.listdir(path) if os.path.isdir(os.path.join(path, elem)) == True if "E_" in elem]
            CVval=[elem.split('_')[-1].strip("/") for elem in glob.glob(path + "**/", recursive=True) if "E_" in elem]

            #sorting by inscreasing order for "levels" in plotting and for sorting colvar files in colvar_list
            CVval2=sorted([float(i) for i in CVval])

            #get dictionary of CVval values with their indexes (last line inverse keys and values)
            indexing = enumerate(CVval2) 
            dico=dict(list(indexing))
            dico={v: k for k, v in dico.items()}

            #reordering values in list according to CVval2 order
            colvar_listS = [None] * len(colvar_list)
            for i in colvar_list:
                colvar_listS[dico[float(re.search(r"(?<=E_).*(?=/)", i).group())]] = i
            colvar_list=colvar_listS

            #cut colvar_list in several lists of 20 elem to lighten the plot generation and avoid crash (for std desktop ressources)
            clcut=[colvar_list[i:i + 20] for i in range(0, len(colvar_list), 20)]

            #function parameters#######################CHANGE CODE TO GET PARAM FROM COLVAR##########################
            a=0.5
            deq=0.05
            x = np.arange(0, 0.2, 0.005)
            y = np.arange(0, 0.2, 0.005)
            X, Y = np.meshgrid(x, y)
            Z = ((deq/X)**a)-((deq/Y)**a)
            #function parameters#######################CHANGE CODE TO GET PARAM FROM COLVAR##########################

            #data objects
            RMSDeq = []
            RMSDax = []
            phi = []
            psi = []

            for nbl, k in enumerate(clcut, 0):
                plt.figure(1)
                plt.figure(figsize=(12,5*len(k))) 
                for pltn, i in zip(range(1, len(k)*2+1, 2), k):

                    RMSDeq = rdata(i)[0]
                    RMSDax = rdata(i)[1]
                    phi = rdata(i)[2]
                    psi = rdata(i)[3]

                    #get param of restraints from the plumed input files, todo: build a dic {colvar_file : CVval} instead of the wierd trick in the following line (should change also in l.87-93)
                    with open(re.search(r".*E_([0-9]|-|\.)*", i).group()+"/plumed_"+format(CVval2[k.index(i)+nbl*20], '.2f')+".dat", "r") as plume:
                        for line in plume:
                            if "KAPPA" in line:
                                param = re.search(r"AT.*KAPPA[^ ]*", line)
                    plume.close()

                    plt.subplot(len(k), 2, pltn)
                    levels=[-1.0, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0]
                    if CVval2[k.index(i)+nbl*20] in levels:
                        levels.remove(CVval2[k.index(i)+nbl*20])
                        CS = plt.contour(X, Y, Z, levels=levels)
                        CS2 = plt.contour(CS, levels=[CVval2[k.index(i)+nbl*20]], colors='red', linestyles='dashed')
                    else:
                        CS = plt.contour(X, Y, Z, levels=levels)
                        CS2 = plt.contour(CS, levels=[CVval2[k.index(i)+nbl*20]], colors='red', linestyles='dashed')
                    plt.clabel(CS, inline=1, fontsize=10)
                    plt.clabel(CS2, inline=1)
                    plt.xlabel('RMSDeq')
                    plt.ylabel('RMSDax')
                    plt.scatter(RMSDeq, RMSDax, c=list(range(1, len(RMSDax)+1, 1)), cmap=plt.cm.get_cmap('rainbow'), s=10)
                    plt.title(param.group())

                    plt.subplot(len(k), 2, pltn+1)
                    plt.xlim(-3.2, 3.2)
                    plt.ylim(-3.2, 3.2)
                    plt.title('phi & psi dihedral angles')
                    plt.xlabel('phi')
                    plt.ylabel('psi')
                    plt.scatter(phi, psi, c=list(range(1, len(RMSDax)+1, 1)), cmap=plt.cm.get_cmap('rainbow'), s=10)


                plt.tight_layout()
                plt.savefig('plot'+str(nbl)+'.jpeg')
                plt.close()

        else:
            print("you should specify the colvar_files with -f or turn umb plot mode on with -umb")

# If called from the command line... 
if __name__ == "__main__":
    main()
