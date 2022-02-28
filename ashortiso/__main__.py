import argparse
import logging as log
from ashortiso.version import __version__
from iqtools import *
import numpy as np
from datetime import datetime

class IsomerIdentification():
    def __init__(self, filename,  inyection_time, time, fcen, fspan, lframes):
        self.fcen=fcen
        self.fspan=fspan
        self.lframes=lframes
        self.itime=inyection_time
        self.tdeath=time
        self.filename=filename
        self.method1()
    
    def _create_spectrogram(self, skip, method=None):
        nframes=int(self.tdeath*self.iq.fs/self.lframes)
        sframes=int(skip*self.iq.fs/self.lframes)
        self.iq.read(nframes=nframes, lframes=self.lframes, sframes=sframes)
        iq.method='mtm' #'fft', 'mtm', 'welch'
        if method: self.iq.method = method
        return self.iq.get_spectrogram(nframes, self.lframes) #f=x[t,p], t=y[p,f], p=z[t,f]
    
    def get_energy_content(self, skip):
        xx,yy,zz= self._create_spectrogram(skip)
        nxx,nyy,nzz=get_cut_spectrogram(xx, yy, zz, xcen=xcen, xspan=xspan)
        axx, ayy, azz = get_averaged_spectrogram(nxx, nyy, nzz, len(nxx[:,0]))
        return np.average(azz[0,:])

    def method1(self):
        self.iq = get_iq_object(self.filename)
        self.iq.read_samples(1)
        energy_isomer= self.get_energy_content(self.itime)
        energy_background=self.get_energy_content(self.itime+2)
        return IsomerIdentification.isomer_or_not(energy_isomer, energy_background)
    
    def method2(self):
        self.iq = get_iq_object(self.filename)
        xx,yy,zz=self._create_spectrogram(0)
        ycen=(self.itime+(self.itime+self.tdeath))/2
        yspan=ycen-self.itime
        nxx,nyy,nzz=get_cut_spectrogram(xx,yy,zz,xcen=xcen, xspan=xspan, ycen=ycen, yspan=yspan)
        bxx,byy,bzz=get_cut_spectrogram(xx,yy,zz,xcen=xcen, xspan=xspan, ycen=ycen, yspan=yspan, invert=True)
        abxx,abyy,abzz=get_averaged_spectrogram(bxx,byy,bzz, len(bxx[:,0]))
        bg=np.average(abzz[0,:])
        nzz=np.substract(nzz,bg)
        
    def frame_by_frame(xx,yy,zz, method='method1'):
        for i,_ in enumerate(yy[:,0]):#time
            xx[i,:] zz[i,:]
                
        
    
        
    @staticmethod
    def isomer_or_not(energy_isomer, energy_background):
        deltaE=np.abs(energy_isomer-energy_background)
        if deltaE > 2*np.sqrt(energy_background): return True
        else: return False

def main():
    scriptname = 'AShortIso'
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='+', help='Name of the input file.')
    parser.add_argument('-t', '--itime', type=float, help='Injection time.')
    parser.add_argument('-td', '--tdeath', type=float, help='Time in which the isomer will have completly dissapeared after the injection time..')    
    parser.add_argument('-fc', '--fcen', type=float,  help='Frecuency center of the supposed isomer.')
    parser.add_argument('-fs', '--fspan', type=float, , help='Frecuency span around fcenter.')
    parser.add_argument('-lf', '--lframes', type=int, default='512' , help='Number of frecuency bins.')
    
    parser.add_argument('-v', '--verbose',
			help='Increase output verbosity', action='store_true')
    parser.add_argument('-out', '--outdir', type=str, default='.',
                                                help='Output directory.')
    
    args = parser.parse_args()

    print(f'Running {scriptname} V{__version__}')
    if args.verbose: log.basicConfig(level=log.DEBUG)
    if args.outdir: outfilepath = os.path.join(args.outdir, '')
    
    # here we go:                                                                                       
    log.info(f'File {args.filename} passed for processing the information of {args.refisotope}+{args.re\
fcharge}.')
    
    date_time=datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    outname=f'{outfilepath}{date_time}-isomers_file.txt'
    with open(f'{outname}','a') as wf:
        files_with_isomers=0
        files_without_isomers=0
        for file in args.filename:
            isomer=IsomerIdentification(file, args.itime, args.tdeath, args.fcen, args.fspan, args.lframes)
            if isomer:
                np.savetxt(wf,f'\x1B[31;1m{file}\x1B[0m', newline='\n') 
                files_with_isomers +=files_with_isomers
            else:
                np.savetxt(wf,f'{file}', newline='\n') 
                files_without_isomers +=files_without_isomers
        total_files_analysed=files_with_isomers+files_without_isomers
        print(f'It has been analysed {total_files_analysed} with isomers present in {files_with_isomers} files and {files_without_isomers} files without evidence of isomers.')
        
if __name__ == '__main__':
    main()
