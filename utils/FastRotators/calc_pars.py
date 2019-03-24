from __future__ import print_function
import pyfits
import argparse
import os
import spfr
import glob

parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-model_path',default='../../data/COELHO2014/')
parser.add_argument('-npools', default=1)
parser.add_argument('-fixG', default=-1)

args   = parser.parse_args()
dirin  = args.directorio
model_path  = args.model_path
npools = int(args.npools)
fixG = float(args.fixG)

try:
    sc = pyfits.getdata(dirin)
    fits = [dirin]
    dirin = dirin.split('/')
    dirin = dirin[:-1]
    dt = ''
    for fn in dirin:
        dt += fn+'/'
    dirin = dt
except:
    fits = glob.glob(dirin + '/*fits')

os.system('mkdir '+dirin+'/FR')
print('\n\n')
for fit in fits:
    print('Estimation of atmospheric parameters for ' + fit + '...')

    sc = pyfits.getdata(fit)
    hd = pyfits.getheader(fit)
    teff, logg, feh, vsini = spfr.get_pars_fr(sc[0],sc[5],model_patht=model_path,npools=npools,fixG=fixG)
    print('\tTeff=', teff, 'K')
    print('\tlog(g)=', logg, 'dex')
    print('\t[Fe/H]=', feh, 'dex')
    print('\tvsini=', vsini, 'km/s\n')


    line_out = "%-15s %18.8f %5.2f %5.2f %5.1f %4.2f %s\n"%(hd['HIERARCH TARGET NAME'], 2400000.5 + float(hd['HIERARCH MBJD']), teff, logg, feh, vsini, fit)

    i    = 0
    isin = False
    if os.access(dirin+'/FR/pars_FR.txt',os.F_OK):
        f = open(dirin+'/FR/pars_FR.txt','r')
        lines = f.readlines()
        for line in lines:
            cos = line.split()
            if cos[1] == line_out.split()[1]:
                lines[i] = line_out
                isin = True
                break
            i+=1
        if not isin:
            lines.append(line_out)
    else:
        lines = [line_out]

    f = open(dirin+'/FR/pars_FR.txt','w')
    for line in lines:
        f.write(line)
    f.close()
