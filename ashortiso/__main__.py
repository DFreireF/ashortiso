import argparse
import logging as log
#from ashortiso.version import __version__
from iqtools import *
import numpy as np
from datetime import datetime
from lmfit import *

class IsomerIdentification():
    def __init__(self, filename,  inyection_time, time, fcen, fspan, lframes):
        self.fcen=fcen
        self.fspan=fspan
        self.itime=inyection_time
        self.tdeath=time
        self.filename=filename
        self.lframes=lframes
    
    def _create_spectrogram(self, skip, lframes, method=None):
        iq = get_iq_object(self.filename)
        iq.read_samples(1)
        nframes=int(self.tdeath*iq.fs/lframes)
        sframes=int(skip*iq.fs/lframes)
        iq.read(nframes=nframes, lframes=self.lframes, sframes=sframes)
        iq.method='mtm' #'fft', 'mtm', 'welch'
        if method: iq.method = method
        return iq.get_spectrogram(nframes, self.lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
    
    @staticmethod
    def fit_gaussian(x,y,amp,cen,wid):
        gmod = Model(IsomerIdentification.gaussian)
        result = gmod.fit(y, x=x, amp=6e-07, cen=0, wid=2e2)
        nxcen=result.params['cen'].value
        return nxcen
    
    def get_energy_content(self, skip, fcen, fspan):
        xx,yy,zz= self._create_spectrogram(skip, self.lframes)
        nxx,nyy,nzz=get_cut_spectrogram(xx, yy, zz, xcen=fcen, xspan=fspan)
        axx, ayy, azz = get_averaged_spectrogram(nxx, nyy, nzz, len(nxx[:,0]))
        return axx, ayy, azz

    
    def method_1(self):
        axx,ayy,azz= self.get_energy_content(self.itime, self.fcen, self.fspan)
        axxb,ayy,azzb=self.get_energy_content(self.itime, -1e4, self.fspan)
        isomer=IsomerIdentification.isomer_or_not(np.average(azz[0,:]), np.average(azzb[0,:]))
        return isomer

    def method2(self):
        xx,yy,zz=self._create_spectrogram(0)
        ycen=(self.itime+(self.itime+self.tdeath))/2
        yspan=ycen-self.itime
        nxx,nyy,nzz=get_cut_spectrogram(xx,yy,zz,xcen=xcen, xspan=xspan, ycen=ycen, yspan=yspan)
        bxx,byy,bzz=get_cut_spectrogram(xx,yy,zz,xcen=xcen, xspan=xspan, ycen=ycen, yspan=yspan, invert=True)
        abxx,abyy,abzz=get_averaged_spectrogram(bxx,byy,bzz, len(bxx[:,0]))
        bg=np.average(abzz[0,:])
        nzz=np.substract(nzz,bg)
        
    @staticmethod
    def isomer_or_not(isomerE, backgroundE):
        deltaE=np.abs(isomerE-backgroundE)
        if deltaE > backgroundE: return True
        else: return False
        
    @staticmethod
    def gaussian(x,amp,cen,wid):
        """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))

def main():
    scriptname = 'AShortIso'
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='+', help='Name of the input file.')
    parser.add_argument('-t', '--itime', type=float, help='Injection time.')
    parser.add_argument('-td', '--tdeath', type=float, help='Time in which the isomer will have completly dissapeared after the injection time..')    
    parser.add_argument('-fc', '--fcen', type=float,  help='Frecuency center of the supposed isomer.')
    parser.add_argument('-fs', '--fspan', type=float, help='Frecuency span around fcenter.')
    parser.add_argument('-lf', '--lframes', type=int, default='512' , help='Number of frecuency bins.')
    
    parser.add_argument('-v', '--verbose',
			help='Increase output verbosity', action='store_true')
    parser.add_argument('-out', '--outdir', type=str, default='.',
                                                help='Output directory.')
    
    args = parser.parse_args()

    print(f'Running {scriptname}') #V{__version__}')
    if args.verbose: log.basicConfig(level=log.DEBUG)
    if args.outdir: outfilepath = os.path.join(args.outdir, '')
    # here we go:
    log.info(f'File {args.filename} passed for processing.')
    
    date_time=datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    outname=f'{outfilepath}{date_time}-isomers_file.txt'
    fwith_iso=0
    fwithout_iso=0
    with open(f'{outname}','a') as wf:
        if ('txt') in args.filename[0]:
            filename_list=read_masterfile(args.filename[0])
            for file in filename_list: fwith_iso, fwithout_iso=files_with_isomers_or_not(file, args.itime, args.tdeath, args.fcen, args.fspan, args.lframes, fwith_iso, fwithout_iso, wf)
        else:
            for file in args.filename: fwith_iso, fwithout_iso=files_with_isomers_or_not(file, args.itime, args.tdeath, args.fcen, args.fspan, args.lframes, fwith_iso, fwithout_iso, wf)
        print_output(fwith_iso, fwithout_iso)
            
def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]

def files_with_isomers_or_not(file, inyec_time, timetostudy, fcen, fspan, binning, fwithiso, fwithoutiso, wf):
    iso=IsomerIdentification(file, args.itime, args.tdeath, args.fcen, args.fspan, args.lframes)
    isomer=iso.method_1()
    if isomer:
        np.savetxt(wf,['##ISO##'+file+'##ISO##'], newline='\n', fmt='%s')
        fwithiso = fwithiso+1
    else:
        np.savetxt(wf, [file], newline='\n', fmt='%s')
        fwithoutiso = fwithoutiso+1
    return fwithiso, fwithoutiso

def print_output(fwith,fwithout):
    total_files_analysed=fwith+fwithout
    print(f'It has been analysed {total_files_analysed} with isomers present in {files_with_isomers} files and {files_without_isomers} files without evidence of isomers.')

if __name__ == '__main__':
    main()
