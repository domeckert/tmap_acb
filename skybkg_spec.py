import sys, getopt
import os
from astropy.io import fits
import urllib.request, urllib.error, urllib.parse
import numpy as np

m1rsp = 'm1.rsp'

def main(argv):
    try:
        opts, args = getopt.getopt(argv, ":p:o:a:b:u:",
                                   ["cxbparfile=", "ofile=", "abund=", "bands=", "uhtp="])
    except getopt.GetoptError:
        print('skybkg_spec.py -p cxbparfile -o outputfile --abund="angr --bands=bands --uhtp=y/n"')
        sys.exit(2)
    if (len(opts) < 2 or len(opts) > 6):
        print('skybkg_spec.py -p cxbparfile -o outputfile --abund="angr" --abund="angr --bands=bands --uhtp=y/n')
        sys.exit(2)
    abund = 'aspl'
    bands = [0.7, 1.2]
    uhtp = 0
    pars = None
    outfile = None
    for opt, arg in opts:
        if opt == '-h':
            print('skybkg_spec.py -p cxbparfile -o outputfile --abund="angr" --abund="angr --bands=bands --uhtp=y/n')
            sys.exit()
        elif opt in ("-o", "--ofile"):
            outfile = arg
        elif opt in ("-p", "--cxbparfile"):
            pars = arg
        elif opt in ("-a", "--abund"):
            abund = arg
        elif opt in ("-b", "--bands"):
            bands = np.loadtxt(arg)
        elif opt in ("-u", "--uhtp"):
            tuh = arg
            if (tuh == 'y' or tuh == 'Y' or tuh == 'yes'):
                uhtp = 1
    if 'HEADAS' not in os.environ:
        print('Error: HEASOFT environment not set, please set it first')
        sys.exit(2)
    # Read parameter file
    if pars is None:
        print('Error: missing mandatory parameter -p')
        sys.exit(2)
    if outfile is None:
        print('Error: missing mandatory parameter -o')
        sys.exit(2)

    fpars = open(pars)
    lp = fpars.readlines()
    fpars.close()

    nh = None
    z = None
    for i in range(len(lp)):
        if 'NH' in lp[i] or 'nh' in lp[i]:
            nh = float(lp[i].split()[1])
        if 'REDSHIFT' in lp[i] or 'z' in lp[i]:
            z = float(lp[i].split()[1])

    if nh is None:
        print('Error: Missing parameter NH')
        sys.exit(2)
    if z is None:
        print('Error: Missing parameter REDSHIFT')
        sys.exit(2)

    lhb, ght, ghn, cxb, uht, uhn = None, None, None, None, None, None
    for i in range(len(lp)):
        if 'lhb' in lp[i]:
            lhb = float(lp[i].split()[1])
            elhb = float(lp[i].split()[2])
        if 'ght' in lp[i]:
            ght = float(lp[i].split()[1])
            eght = float(lp[i].split()[2])
        if 'ghn' in lp[i]:
            ghn = float(lp[i].split()[1])
            eghn = float(lp[i].split()[2])
        if 'cxb' in lp[i]:
            cxb = float(lp[i].split()[1])
            ecxb = float(lp[i].split()[2])
        if 'uht' in lp[i]:
            uht = float(lp[i].split()[1])
            euht = float(lp[i].split()[2])
        if 'uhn' in lp[i]:
            uhn = float(lp[i].split()[1])
            euhn = float(lp[i].split()[2])

    if lhb is None:
        print('Error: missing parameter lhb in input file')
        sys.exit(2)
    if ght is None:
        print('Error: missing parameter ght in input file')
        sys.exit(2)
    if ghn is None:
        print('Error: missing parameter ghn in input file')
        sys.exit(2)
    if cxb is None:
        print('Error: missing parameter cxb in input file')
        sys.exit(2)
    if uhtp:
        if uht is None:
            print('Error: missing parameter uht in input file')
            sys.exit(2)
        if uhn is None:
            print('Error: missing parameter uhn in input file')
            sys.exit(2)

    # Open output XCM file
    fxcm = open(outfile+'.xcm', 'w')
    fxcm.write('statistic cstat\n')
    fxcm.write('\n')
    fxcm.write('cosmo 70 0 0.7\n')
    fxcm.write('abund %s\n' % (abund))
    fxcm.write('\n')

    # Write sky model
    fitbkg = False
    if (fitbkg):
        if not uhtp and not pspc:
            fxcm.write('model constant(apec + phabs(apec + powerlaw))\n')
        elif uhtp and pspc:
            fxcm.write('model constant(apec + phabs(apec + powerlaw + apec + apec))\n')
        else:
            fxcm.write('model constant(apec + phabs(apec + powerlaw + apec))\n')
    else:
        if not uhtp:
            fxcm.write('model constant(apec + phabs(apec + powerlaw + apec))\n')
        else:
            fxcm.write('model constant(apec + phabs(apec + powerlaw + apec + apec))\n')

    # Write sky model
    fxcm.write('1.0      -0.01          0          0      1e+10      1e+10\n')
    fxcm.write('0.11         -1      0.008      0.008         64         64\n')
    fxcm.write('1     -0.001          0          0          5          5\n')
    fxcm.write('0      -0.01          0          0         10         10\n')
    fxcm.write('%g     -0.01         %g         %g      %g      %g\n' % (lhb, lhb - elhb, lhb + elhb, lhb - elhb, lhb + elhb))
    fxcm.write('%g     -0.01          0          0     100000      1e+06\n' % (nh))
    fxcm.write('%g      -0.01         %g         %g      %g      %g\n' % (ght, ght - eght, ght + eght, ght - eght, ght + eght))
    fxcm.write('1     -0.001          0          0          5          5\n')
    fxcm.write('0      -0.01          0          0         10         10\n')
    fxcm.write('%g      -0.01         %g         %g      %g      %g\n' % (ghn, ghn - eghn, ghn + eghn, ghn - eghn, ghn + eghn))
    fxcm.write('1.46         -1         -3         -2          9         10\n')
    fxcm.write('%g       -0.01         %g         %g      %g      %g\n' % (cxb, cxb - ecxb, cxb + ecxb, cxb - ecxb, cxb + ecxb))
    fxcm.write('5.0      -0.01      0.2      0.2         25         25\n')
    fxcm.write('0.3         0.01          0.01          0.01         1.5          1.5\n')
    fxcm.write('%g      -0.01          0          0         10         10\n' % (z))
    fxcm.write('0      -1e-07          0          0      10      10\n')
    if uhtp:
        fxcm.write('%g       -0.01         %g         %g      %g      %g\n' % (uht, uht - euht, uht + euht, uht - euht, uht + euht))
        fxcm.write('1     -0.001          0          0          5          5\n')
        fxcm.write('0      -0.01          0          0         10         10\n')
        fxcm.write('%g       -0.01         %g         %g      %g      %g\n' % (uhn, uhn - euhn, uhn + euhn, uhn - euhn, uhn + euhn))

    if os.path.exists('m1.fak'):
        os.system('rm m1.fak')

    fxcm.write('fakeit none\n')
    fxcm.write('%s\n' % (m1rsp))
    fxcm.write('\n')
    fxcm.write('n\n')
    fxcm.write('\n')
    fxcm.write('m1.fak\n')
    fxcm.write('1000, 1, 1000\n')

    print('Hello')
    lowthresh=bands[:,0]
    highthresh=bands[:,1]
    print(lowthresh)
    print(highthresh)
    nband = len(lowthresh)
    ratefile = outfile + '_rate.txt'
    for nb in range(nband):
        fxcm.write('ign **-%1.3lf %1.3lf-**\n' % (lowthresh[nb], highthresh[nb]))
        fxcm.write('tclout rate 1\n')
        fxcm.write('echo "Band %g %g: $xspec_tclout" >> %s\n' % (lowthresh[nb], highthresh[nb], ratefile))
        fxcm.write('not 0.3-12.0\n')
    	
    fxcm.write('quit\n')
    fxcm.close()

    os.system('xspec < %s' % (outfile+'.xcm'))

if __name__ == "__main__":
    main(sys.argv[1:])
