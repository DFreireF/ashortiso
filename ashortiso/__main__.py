import argparse
import logging as log
from ashortiso.version import __version__
from iqtools import *
import numpy as np

class IsomerIdentification():
    def __init__(self, filename,  inyection_time, time, fcen, fspan, lframes):
        self.fcen=fcen
        self.fspan=fspan
        self.lframes=lframes
        self.itime=inyection_time
        self.tdeath=time
        self.filename=filename
        
    def get_energy_content(self, skip):
        xx,yy,zz= self._create_spectrogram(skip)
        nxx,nyy,nzz=get_cut_spectrogram(xx, yy, zz, xcen=xcen, xspan=xspan)
        axx, ayy, azz = get_averaged_spectrogram(nxx, nyy, nzz, len(nxx[:,0]))
        return np.average(azz[0,:])
    
    def _create_spectrogram(self, skip, method=None):
        nframes=int(self.tdeath*self.iq.fs/self.lframes)
        sframes=int(skip*self.iq.fs/self.lframes)
        self.iq.read(nframes=nframes, lframes=self.lframes, sframes=sframes)
        iq.method='mtm' #'fft', 'mtm', 'welch'
        if method: self.iq.method = method
        return self.iq.get_spectrogram(nframes, self.lframes) #f=x[t,p], t=y[p,f], p=z[t,f]

    def method1(self):
        self.iq = get_iq_object(self.filename)
        self.iq.read_samples(1)
        energy_isomer= self.get_energy_content(self.itime)
        energy_background=self.get_energy_content(self.itime+2)
        return np.abs(energy_isomer-energy_background)
    
    @staticmethod
    def isomer_or_not(energy_isomer, energy_background):
        deltaE=np.abs(energy_isomer-energy_background)
        if deltaE > 1.5*energy_background: return True
        else: return False

def main():
    scriptname = 'ashortiso'
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='+', help='Name of the input file.')
    parser.add_argument('-t', '--itime', type=float, help='Injection time.')
    parser.add_argument('-td', '--tdeath', type=float, help='Time in which the isomer will have completly dissapeared after the injection time..')    
    parser.add_argument('-fc', '--fcen', type=float,  help='Frecuency center of the supposed isomer.')
    parser.add_argument('-fs', '--fspan', type=float, , help='Frecuency span around fcenter.')
    parser.add_argument('-lf', '--lframes', type=int, default='512' , help='Number of frecuency bins.')
    
    parser.add_argument('-v', '--verbose',
			help='Increase output verbosity', action='store_true')

    args = parser.parse_args()

    print(f'Running {scriptname} V{__version__}')
    if args.verbose: log.basicConfig(level=log.DEBUG)

    # here we go:                                                                                       
    log.info(f'File {args.filename} passed for processing the information of {args.refisotope}+{args.re\
fcharge}.')
    
    for file in args.filename:
        isomer=IsomerIdentification(file, args.itime, args.tdeath, args.fcen, args.fspan, args.lframes)
        files_with_isomers=0
        files_without_isomers=0
        if isomer:
            save namefile '\x1B[31;1m' '\x1B[0m'
            files_with_isomers +=files_with_isomers
        else:
            save namefile
            files_without_isomers +=files_with_isomers
            
        total_files_analysed=files_with_isomers+files_without_isomers
        print(f'It has been analysed {total_files_analysed} with isomers present in {files_with_isomers} files and {files_without_isomers} files without evidence of isomers.')
        

if __name__ == '__main__':
    main()
