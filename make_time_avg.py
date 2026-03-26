import numpy as np
import os


def find_nlast(dirr,nfirst=0):
    i = nfirst
    while os.path.exists(dirr+f'/gasdens{i}.dat'):
        i+=1
    return i

def get_domain(dirr):
    X = np.loadtxt(dirr+'/domain_x.dat')
    nx = len(X)-1
    Y = np.loadtxt(dirr+'/domain_y.dat')[3:-3]
    ny = len(Y)-1
    Z = np.loadtxt(dirr+'/domain_z.dat')[3:-3]
    nz = len(Z)-1
    shape = (nz,ny,nx)
    return X,Y,Z,shape

def get_state_avg(dirr,name,nlast=None,nfirst=0):
    X,Y,Z,shape = get_domain(dirr)
    avgarr = np.zeros(shape)
    if nlast is None:
        nlast = find_nlast(dirr,nfirst)
    N = nlast-nfirst
    for i in range(nfirst,nlast):
        arr = np.fromfile(dirr+f'/{name}{i}.dat').reshape(shape)
        avgarr += arr
    return avgarr/N

def write_planet_avg_file(dirr1,dirr2,nlast=None,nfirst=0):
    if nlast is None:
        nlast = find_nlast(dirr1,nfirst)
    avgx=0
    avgy=0
    avgz=0
    avgvx=0
    avgvy=0
    avgvz=0
    avgmass=0
    avgtime=0
    avgomframe=0
    N = 0
    with open(dirr1+'/planet0.dat') as f:
        for line in f:
            nout = int(line.split()[0])
            if nout < nfirst or nout >= nlast:
                continue
            nout,x,y,z,vx,vy,vz,mass,time,omframe = map(float,line.split())
            avgx+=x
            avgy+=y
            avgz+=z
            avgvx+=vx
            avgvy+=vy
            avgvz+=vz
            avgmass+=mass
            avgtime+=time
            avgomframe+=omframe
            N += 1
    avgx/=N
    avgy/=N
    avgz/=N
    avgvx/=N
    avgvy/=N
    avgvz/=N
    avgmass/=N
    avgtime/=N
    avgomframe/=N
    with open(dirr2+'/planet0.dat','w+') as f:
        f.write(f'avg\t{avgx}\t{avgy}\t{avgz}\t{avgvx}\t{avgvy}\t{avgvz}\t{avgmass}\t{avgtime}\t{avgomframe}\n')
    return 0

def copy_domain_files(dirr1,dirr2):
    # copy domain files from dirr1 to dirr2
    import shutil
    shutil.copy(dirr1+'/domain_x.dat',dirr2)
    shutil.copy(dirr1+'/domain_y.dat',dirr2)
    shutil.copy(dirr1+'/domain_z.dat',dirr2)
    return 0

def main():
    fulldirr = '/project2/fciesla/ericvc/fargo/outputs/alpha4_mplan300'
    outdirr = '/project2/fciesla/ericvc/fargo/outputs/alpha4_mplan300_115-125'
    if not os.path.exists(outdirr):
        print('mkdir ',outdirr)
        os.makedirs(outdirr)
    X,Y,Z,shape = get_domain(fulldirr)
    #nlast = find_nlast(fulldirr)
    nlast = 125
    print(f'{nlast = }')
    nfirst = 115

    # copy domain files
    copy_domain_files(fulldirr,outdirr)

    # get time average state variables
    for state in ['gasdens','gasenergy','gasvx','gasvy','gasvz']:
        print('working on ',state)
        avgarr = get_state_avg(fulldirr,state,nlast,nfirst=nfirst)
        #np.save(outdirr+f'/{state}avg.dat',avgarr)
        avgarr.tofile(outdirr+f'/{state}avg.dat')

    # copy variables file
    import shutil
    shutil.copy(fulldirr+'/variables.par',outdirr)

    # copy summary file
    shutil.copy(fulldirr+'/summary0.dat',outdirr+'/summaryavg.dat')

    # copy planet data
    write_planet_avg_file(fulldirr,outdirr,nlast)
    
    

if __name__ == '__main__':
    main()

