import sys, getopt
import os
from astropy.io import fits
import urllib.request, urllib.error, urllib.parse
import numpy as np

rsp = 'm1.rsp'

def main(argv):
    try:
        opts, args = getopt.getopt(argv, ":z:n:o:a:b:w:",
                                   ["redshift=", "nh=", "ofile=", "abund=", "bands=", "withabund="])
    except getopt.GetoptError:
        print('templates.py -z redshift -n nh -o outputfile --abund="angr" --bands=bands --withabund=y/n"')
        sys.exit(2)
    if (len(opts) < 2 or len(opts) > 6):
        print('templates.py -z redshift -n nh -o outputfile --abund="angr" --bands=bands --withabund=y/n')
        sys.exit(2)
    abund = 'aspl'
    bands = [0.7, 1.2]
    nh = None
    z = None
    pars = None
    outfile = None
    withabund = False
    for opt, arg in opts:
        if opt == '-h':
            print('templates.py -z redshift -n nh -o outputfile --abund="angr" --bands=bands --withabund=y/n')
            sys.exit()
        elif opt in ("-o", "--ofile"):
            outfile = arg
        elif opt in ("-z", "--redshift"):
            z = float(arg)
        elif opt in ("-n", "--nh"):
            nh = float(arg)
        elif opt in ("-a", "--abund"):
            abund = arg
        elif opt in ("-b", "--bands"):
            bands = np.loadtxt(arg)
        elif opt in ("-w", "--withabund"):
            if arg=='y':
                withabund = True

    if 'HEADAS' not in os.environ:
        print('Error: HEASOFT environment not set, please set it first')
        sys.exit(2)
    # Read parameter file
    if z is None:
        print('Error: missing mandatory parameter -z')
        sys.exit(2)
    if nh is None:
        print('Error: missing mandatory parameter -n')
        sys.exit(2)
    if outfile is None:
        print('Error: missing mandatory parameter -o')
        sys.exit(2)

    nkt = 50
    nab = 20

    kt_grid = np.logspace(np.log10(0.2), np.log10(25), nkt)
    if withabund:
        ab_grid = np.linspace(0, 1.5, nab)
    else:
        ab_grid = np.array([0.3])
        
    #kt_grid = [5, 2]
    #ab_grid = [0.3, 0.4]

    lowthresh=bands[:,0]
    highthresh=bands[:,1]
    print(lowthresh)
    print(highthresh)

    nband = len(lowthresh)

    fout = open(outfile, 'w')
    fout.write('# kt ab')
    for band in range(nband):
        tst = str(lowthresh[band])+'_'+str(highthresh[band])
        fout.write(' %s' % (tst))

    fout.write('\n')

    rsproot = rsp.split('.')[0]

    for ikt, kt in enumerate(kt_grid):

        for iab, ab in enumerate(ab_grid):

            simout = rsproot+'.fak'

            # Open output XCM file
            fxcm = open('commands.xcm', 'w')
            fxcm.write('statistic cstat\n')
            fxcm.write('\n')
            fxcm.write('cosmo 70 0 0.7\n')
            fxcm.write('abund %s\n' % (abund))
            fxcm.write('\n')
            fxcm.write('model phabs(apec)\n')
            fxcm.write('%g     -0.01          0          0     100000      1e+06\n' % (nh))
            fxcm.write('%g      -0.01      0.15      0.15         26         26\n' % (kt))
            fxcm.write('%g      -0.01      0      0         2         2\n' % (ab))
            fxcm.write('%g      -0.01      0      0         10         10\n' % (z))
            fxcm.write('1     -0.001          0          0          5          5\n')

            if ikt==0 and iab==0:
                if os.path.exists(simout):
                    os.remove(simout)

                fxcm.write('fakeit none\n')
                fxcm.write('%s\n' % (rsp))
                fxcm.write('\n')
                fxcm.write('n\n')
                fxcm.write('\n')
                fxcm.write('%s\n' % (simout))
                fxcm.write('\n')

            else:
                fxcm.write('dat 1:1 %s\n' % (simout))

            ratefile = outfile + '_rate.txt'
            if os.path.exists(ratefile):
                os.remove(ratefile)

            for nb in range(nband):
                fxcm.write('ign **-%1.3lf %1.3lf-**\n' % (lowthresh[nb], highthresh[nb]))
                fxcm.write('tclout rate 1\n')
                fxcm.write('echo "Band %g %g: $xspec_tclout" >> %s\n' % (lowthresh[nb], highthresh[nb], ratefile))
                fxcm.write('not 0.3-12.0\n')

            fxcm.write('quit\n')
            fxcm.close()

            os.system('xspec < commands.xcm')

            if withabund:
                fout.write('%g  %g  '% (kt, ab))
            else:
                fout.write('%g '% (kt))

            fr = open(ratefile)
            lr = fr.readlines()
            fr.close()

            for nb in range(nband):

                cr = float(lr[nb].split()[5])
                fout.write('%1.4lf ' % (cr))

            fout.write('\n')


if __name__ == "__main__":
    main(sys.argv[1:])
