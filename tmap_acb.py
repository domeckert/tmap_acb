cimport numpy as np
from astropy.io import fits
import sys,getopt
import iminuit
from scipy.interpolate import interp1d
from astropy.wcs import WCS

def extract_skybkg(all_img, all_exp, cx, cy, rmax, pixsize, all_bkg=None):

    nimg = len(all_img)
    shape = all_img[0].shape
    
    y,x = np.indices(shape)
    rads = np.hypot(y-cy, x-cx) * pixsize # arcmin
    
    skybkg, skybkg_err = np.empty(nimg), np.empty(nimg)

    for i in range(nimg):
        
        timg = all_img[i]
        texp = all_exp[i]
        
        sel_reg = np.where(np.logical_and(rads>rmax,texp>0))
        
        nv = len(sel_reg[0])
            
        if all_bkg is not None:
        
            tbkg = all_bkg[i]
            
            nxb = np.sum(tbkg[sel_reg] / texp[sel_reg]) / nv / pixsize ** 2
        else:
        
            nxb = 0.
            
        skybkg[i] = np.sum(timg[sel_reg] / texp[sel_reg]) / nv / pixsize ** 2 - nxb    
    
        skybkg_err[i] = np.sqrt(np.sum(timg[sel_reg] / texp[sel_reg] ** 2)) / nv / pixsize ** 2
        
    return skybkg, skybkg_err
    
def region(exposure, wcs_inp, pixsize, regfile):
    '''
    Filter out regions provided in an input DS9 region file

    :param regfile: Path to region file. Accepted region file formats are fk5 and image.
    :type regfile: str
    '''
    freg = open(regfile)
    lreg = freg.readlines()
    freg.close()
    nsrc = 0
    nreg = len(lreg)
    expo = np.copy(exposure)
    axes = exposure.shape
    y, x = np.indices(axes)
    regtype = None
    for i in range(nreg):
        if 'fk5' in lreg[i]:
            regtype = 'fk5'
        elif 'image' in lreg[i]:
            regtype = 'image'
    if regtype is None:
        print('Error: invalid format')
        return
    for i in range(nreg):
        if 'circle' in lreg[i]:
            vals = lreg[i].split('(')[1].split(')')[0]
            if regtype == 'fk5':
                xsrc = float(vals.split(',')[0])
                ysrc = float(vals.split(',')[1])
                rad = vals.split(',')[2]
                if '"' in rad:
                    rad = float(rad.split('"')[0]) / pixsize / 60.
                elif '\'' in rad:
                    rad = float(rad.split('\'')[0]) / pixsize
                else:
                    rad = float(rad) / pixsize * 60.
                wc = np.array([[xsrc, ysrc]])
                pixcrd = wcs_inp.wcs_world2pix(wc, 1)
                xsrc = pixcrd[0][0] - 1.
                ysrc = pixcrd[0][1] - 1.
            else:
                xsrc = float(vals.split(',')[0])
                ysrc = float(vals.split(',')[1])
                rad = float(vals.split(',')[2])

            # Define box around source to spped up calculation
            boxsize = np.round(rad + 0.5).astype(int)
            intcx = np.round(xsrc).astype(int)
            intcy = np.round(ysrc).astype(int)
            xmin = np.max([intcx-boxsize, 0])
            xmax = np.min([intcx+boxsize + 1, axes[1]])
            ymin = np.max([intcy-boxsize, 0])
            ymax = np.min([intcy+boxsize + 1, axes[0]])
            rbox = np.hypot(x[ymin:ymax,xmin:xmax] - xsrc,y[ymin:ymax,xmin:xmax] - ysrc)
            # Mask source
            src = np.where(rbox < rad)
            expo[ymin:ymax,xmin:xmax][src] = 0.0
            nsrc = nsrc + 1
        elif 'ellipse' in lreg[i]:
            vals = lreg[i].split('(')[1].split(')')[0]
            if regtype == 'fk5':
                xsrc = float(vals.split(',')[0])
                ysrc = float(vals.split(',')[1])
                rad1 = vals.split(',')[2]
                rad2 = vals.split(',')[3]
                angle = float(vals.split(',')[4])
                if '"' in rad1:
                    rad1 = float(rad1.split('"')[0]) / pixsize / 60.
                    rad2 = float(rad2.split('"')[0]) / pixsize / 60.
                elif '\'' in rad1:
                    rad1 = float(rad1.split('\'')[0]) / pixsize
                    rad2 = float(rad2.split('\'')[0]) / pixsize
                else:
                    rad1 = float(rad1) / pixsize * 60.
                    rad2 = float(rad2) / pixsize * 60.
                wc = np.array([[xsrc, ysrc]])
                pixcrd = wcs_inp.wcs_world2pix(wc, 1)
                xsrc = pixcrd[0][0] - 1.
                ysrc = pixcrd[0][1] - 1.
            else:
                xsrc = float(vals.split(',')[0])
                ysrc = float(vals.split(',')[1])
                rad1 = float(vals.split(',')[2])
                rad2 = float(vals.split(',')[3])
                angle = float(vals.split(',')[2])
            ellang = angle * np.pi / 180. + np.pi / 2.
            aoverb = rad1/rad2
            # Define box around source to spped up calculation
            boxsize = np.round(np.max([rad1, rad2]) + 0.5).astype(int)
            intcx = np.round(xsrc).astype(int)
            intcy = np.round(ysrc).astype(int)
            xmin = np.max([intcx-boxsize, 0])
            xmax = np.min([intcx+boxsize + 1, axes[1]])
            ymin = np.max([intcy-boxsize, 0])
            ymax = np.min([intcy+boxsize + 1, axes[0]])
            xtil = np.cos(ellang)*(x[ymin:ymax,xmin:xmax]-xsrc) + np.sin(ellang)*(y[ymin:ymax,xmin:xmax]-ysrc)
            ytil = -np.sin(ellang)*(x[ymin:ymax,xmin:xmax]-xsrc) + np.cos(ellang)*(y[ymin:ymax,xmin:xmax]-ysrc)
            rbox = aoverb * np.hypot(xtil, ytil / aoverb)
            # Mask source
            src = np.where(rbox < rad1)
            expo[ymin:ymax,xmin:xmax][src] = 0.0
            nsrc = nsrc + 1

    print('Excluded %d sources' % (nsrc))
    return expo


def fit_values(all_img, all_exp, fint, cx, cy, ncount, pixsize, skybkg, all_bkg=None):
    
    nphot = 0
    rad_pix = 0.
    
    imgsel = all_img[0]
    
    y, x = np.indices(imgsel.shape)
    rads = np.hypot(y-cy, x-cx)

    while nphot < ncount:

        rad_pix = rad_pix + 0.5

        sel = np.where(np.logical_and(rads<rad_pix, all_exp[0]>0))

        nphot = np.sum(imgsel[sel])

    nv = len(sel[0])
    
    nimg = len(all_img)
    
    counts, bkgcounts, effexp = np.zeros(nimg), np.zeros(nimg), np.zeros(nimg)
    
    area = nv * pixsize ** 2
    
    for i in range(nimg):
    
        counts[i] = np.sum(all_img[i][sel])
        
        if all_bkg is not None:
        
            bkgcounts[i] = np.sum(all_bkg[i][sel])
            
        effexp[i] = np.sum(all_exp[i][sel]) / nv
        
    def calc_cstat(kt, norm):
    
        model = np.empty(nimg)
        
        for i in range(nimg):
        
            tinterp = fint[i](kt)
            
            if tinterp <= 0. or kt<=0. or norm<=0.:
            
                model[i] = skybkg[i] * area * effexp[i] + bkgcounts[i]
                
            else:
            
                model[i] = (norm * tinterp + skybkg[i]) * area * effexp[i] + bkgcounts[i]
                        
        cstat = 0.
        
        for i in range(nimg):
        
            if counts[i]>0:
            
                cstat = cstat + 2. * (model[i] - counts[i] * np.log(model[i]) - counts[i] + counts[i] * np.log(counts[i])) # normalized C-statistic
                
            else:
            
                cstat = cstat + 2. * model[i] 

        return cstat

    minuit = iminuit.Minuit(calc_cstat, kt=5., norm=1e-4)

    minuit.errordef = 1

    out = minuit.migrad()

    kt = out.values['kt']
    
    ekt = out.errors['kt']
    
    norm = out.values['norm']
    
    enorm = out.errors['norm']
                
    return kt, ekt, norm, enorm, rad_pix * pixsize

       

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"i:e:t:o:m:c:b:r:i:s:g:",["images=","expmaps=","template=","mask=","outfile=","ncount=","bkgmaps=","rmax=","tmin=","regfile=","skybkg="])
    except getopt.GetoptError:
        print('tmap_acb.py -i imglist -e explist -t template -o outfile --bkgmaps=bkglist --mask=mask.fits --ncount=200 --rmax=3 --tmin=0.5 --regfile=file.reg --skybkg=skybkg.dat')
        sys.exit(2)
    if (len(opts)<1 or len(opts)>12):
        print('tmap_acb.py -i imglist -e explist -t template -o outfile --bkgmaps=bkglist --mask=mask.fits --ncount=200 --rmax=3 --tmin=0.5 --regfile=file.reg --skybkg=skybkg.dat')
        sys.exit(2)


        # Read arguments
    imglist, explist, bkglist = None, None, None
    outfile = None
    template_file = None
    maskfile = None
    ncount = 200
    rmax = 3. # maximum radius in arcmin outside of which we compute the sky background
    tmin = None # temperatures below this value will be cut
    regfile = None
    fit_skybkg = True
    skybkg_file = None
    for opt, arg in opts:
        if opt in ("-i", "--images"):
            imglist = arg
            print('Image file list is', imglist)
        elif opt in ("-e", "--expmaps"):
            explist = arg
            print('Exposure map file list is', explist)
        elif opt in ("-b", "--bkgmaps"):
            bkglist = arg
            print('Background map file list is', bkglist)
        elif opt in ("-t", "--template"):
            template_file = arg
            print('Template file is', template_file)
        elif opt in ("-o", "--outfile"):
            outfile = arg
            print('Output file list is', outfile)
        elif opt in ("-m", "--mask"):
            maskfile = arg
            print('Mask file is', maskfile)
        elif opt in ("-c", "--ncount"):
            ncount = int(arg)
            print('We will use a target of',ncount,'count per bin')
        elif opt in ("-r", "--rmax"):
            rmax = float(arg)
            print('We will determine the sky background spectrum at R>',rmax,'from the image center')
        elif opt in ("-i", "--tmin"):
            tmin = float(arg)
            print('We will cut temperatures below',tmin)
        elif opt in ("-r", "--regfile"):
            regfile = arg
            print('We will filter out sources in file',regfile)
        elif opt in ("-g", "--skybkg"):
            skybkg_file = arg
            print('We will use the sky background count rates provide in file',skybkg_file)           


    if imglist is None:
        print('No image list provided, aborting')
        sys.exit(2)
        
    if explist is None:
        print('No exposure map list provided, aborting')
        sys.exit(2)
            
    if template_file is None:
        print('No template file provided, aborting')
        sys.exit(2)
            
    if outfile is None:
        print('No output file name provided, we will write the result to default file name tmap.fits')
        outfile = 'tmap.fits'
        

    fimg = open(imglist)
    limg = fimg.readlines()
    fimg.close()

    nimg = len(limg)

    fexp = open(explist)
    lexp = fexp.readlines()
    fexp.close()

    if len(lexp)!=nimg:
        print('Error: the number of provided exposure maps is different from the number of input images')
        sys.exit(2)
        
    if bkglist is not None:
        fbkg = open(bkglist)
        lbkg = fbkg.readlines()
        fbkg.close()

        if len(lbkg)!=nimg:
            print('Error: the number of provided background maps is different from the number of input images')
            sys.exit(2)
        
    template = np.loadtxt(template_file)

    if template.shape[1] != nimg + 1:

        print('Error: the number of provided energy band templates is different from the number of input images')

    all_img, all_exp, all_bkg = [], [], []
    
    pixsize, header = None, None
    
    shape = None

    wcs_inp = None
    
    for i in range(nimg):
        
        tf = limg[i].split()[0]
        
        print('Reading image %d: %s' % (i+1, tf))
        fimg = fits.open(tf)
        timg = fimg[0].data
                
        if i == 0:

            shape = timg.shape

            pixsize = fimg[0].header['CDELT2']*60.
            header = fimg[0].header
            wcs_inp = WCS(header, relax=False)

        else:
        
            temp_shape = timg.shape
            if temp_shape != shape:
        
                print('Error: image shape is different from previous images')
                sys.exit(2)
        
        all_img.append(timg)
        fimg.close()

        tf = lexp[i].split()[0]
        print('Reading exposure map %d: %s' %(i+1,tf))        
        fexp = fits.open(tf)
        texp = fexp[0].data
        temp_shape = texp.shape
        if temp_shape != shape:
    
            print('Error: exposure map shape is different from previous images')
            sys.exit(2)

        if regfile is not None:
            expo = region(texp, wcs_inp, pixsize, regfile)
            all_exp.append(expo)

        else:
            all_exp.append(texp)

        fexp.close()
        
        if bkglist is not None:
        
            tf = lbkg[i].split()[0]
            print('Reading background map %d: %s' % (i+1,tf))
            fbkg = fits.open(tf)
            tbkg = fbkg[0].data
            temp_shape = tbkg.shape
            if temp_shape != shape:
        
                print('Error: background map shape is different from previous images')
                sys.exit(2)
        
            all_bkg.append(tbkg)
            fbkg.close()
        
        
    if maskfile is not None:
        fmask = fits.open(maskfile)
        mask = fmask[0].data
        temp_shape = mask.shape
        
        if temp_shape != shape:
            print('Error: mask shape is different from image shape')
            sys.exit(2)
            
        fmask.close()
        
    cx, cy = shape[1]/2., shape[0]/2.

    y, x = np.indices(shape)
    rads = np.hypot(y - cy, x - cx) * pixsize  # arcmin

    if fit_skybkg:
        skybkg, skybkg_err = extract_skybkg(all_img=all_img,
                                            all_exp=all_exp,
                                            cx=cx,
                                            cy=cy,
                                            rmax=rmax,
                                            pixsize=pixsize,
                                            all_bkg=all_bkg)

    else:
        skybkg = np.empty(nimg)
        fsb = open(skybkg_file)
        lsb = fsb.readlines()
        fsb.close()

        for i in range(nimg):
            skybkg[i] = float(lsb[i].split()[3])
                                        
    all_fint = []
    for i in range(1,nimg+1):
    
        fint = interp1d(template[:,0], template[:,i], kind='cubic', fill_value='extrapolate')
        
        all_fint.append(fint)
        
    if tmin is None:
    
        tmin = np.min(template[:,0])
        
    tmap, nmap, ektmap, enmap, radmap = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(shape)

    if maskfile is None:

        mask = np.ones(shape)


    zer = np.where(all_exp[0]==0)

    mask[zer] = 0.

    for j in range(shape[0]):
        
        print('Row ',j)
        
        for i in range(shape[1]):
            
            if mask[j, i] > 0 and rads[j, i]<=rmax:

                tkt, ekt, norm, enorm, rad = fit_values(all_img=all_img,
                                                          all_exp=all_exp,
                                                          fint=all_fint,
                                                          cx=i,
                                                          cy=j,
                                                          ncount=ncount,
                                                          pixsize=pixsize,
                                                          skybkg=skybkg,
                                                          all_bkg=all_bkg)

                tmap[j,i] = tkt

                ektmap[j,i] = ekt

                nmap[j,i] = norm

                enmap[j,i] = enorm

                radmap[j,i] = rad

    lowT = np.where(tmap<tmin)
    tmap[lowT] = 0.
    ektmap[lowT] = 0.
    nmap[lowT] = 0.
    enmap[lowT] = 0.
    radmap[lowT] = 0.
    lowsnr = np.where(nmap < 2. * enmap)
    tmap[lowsnr] = 0.
    ektmap[lowsnr] = 0.
    nmap[lowsnr] = 0.
    enmap[lowsnr] = 0.
    radmap[lowsnr] = 0.

    hdul = fits.HDUList([fits.PrimaryHDU()])

    tmaphdu = fits.ImageHDU(tmap, name='KT', header=header)
    etmaphdu = fits.ImageHDU(ektmap, name='ERR_KT', header=header)
    normhdu = fits.ImageHDU(nmap, name='NORM', header=header)
    enormhdu = fits.ImageHDU(enmap, name='ERR_NORM', header=header)
    radhdu = fits.ImageHDU(radmap, name='RADIUS', header=header)
    phdu = fits.ImageHDU(tmap*np.nan_to_num(np.sqrt(nmap)), name='PRESSURE', header=header)
    ephdu = fits.ImageHDU(ektmap*np.nan_to_num(np.sqrt(nmap)), name='ERR_P', header=header)
    Khdu = fits.ImageHDU(np.nan_to_num(tmap*np.power(nmap,-1./3.)), name='ENTROPY', header=header)
    eKhdu = fits.ImageHDU(np.nan_to_num(ektmap*np.power(nmap,-1./3.)), name='ERR_K', header=header)

    hdul.append(tmaphdu)
    hdul.append(etmaphdu)
    hdul.append(normhdu)
    hdul.append(enormhdu)
    hdul.append(radhdu)
    hdul.append(phdu)
    hdul.append(ephdu)
    hdul.append(Khdu)
    hdul.append(eKhdu)
    
    hdul.writeto(outfile, overwrite=True)

if __name__ == "__main__":
    main(sys.argv[1:])


