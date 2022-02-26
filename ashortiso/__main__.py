import argparse
import logging as log
from ashortiso.version import __version__
from iqtools import *


def _cut_spectrogram(spec_to_cut, f1, f2, y1=None, y2=None):
    xcen=(f1+f2-2*self.iq.center)/2
    xspan=abs(f1-f2)
    nxx, nzz, nyy= get_cut_spectrogram(spect_to_cut[0], spec_to_cut[1], spec_to_cut[2], xcen=xcen, xspan=xspan)
    return (np.stack((nxx, nzz, nyy), axis=1)).reshape((len(nxx),3))
    
def _create_spectrogram(self, howlong_t, where_t, lf, method=None):
    nframes=int(howlong_t*self.iq.fs/lf)
    sframes=int(where_t*self.iq.fs/lf)
    self.iq.read(nframes=nframes, lframes=lf, sframes=sframes)
    iq.method='mtm' #'fft', 'mtm', 'welch'                                                                                                                                                                     
    if method: iq.method = method
    xx, yy, zz = self.iq.get_spectrogram(nframes,lframes)
    return (np.stack((xx, zz, yy), axis=1)).reshape((len(xx),3))

def method1(filename, lframes, inject_t):
    self.iq = get_iq_object(filename)
    self.iq.read_samples(1)
    time=inject_t-2
    skip_time=0
    bg1 = self._create_spectrogram(time, skip_time, lframes)
    nbg1= self._cut_spectrogram(bg1, f1, f2)

    time=0.100 #100ms of data                                                                                                                                                                                  
    skip_time=inject_t
    interest_data = self._create_spectrogram(time, skip_time, lframes)
    idata= self._cut_spectrogram(interest_data, f1, f2)
    
    time=inject_t-2
    skip_time=inject_t+1
    bg2 = ImportData._create_spectrogram(time, skip_time, lframes)
    nbg2= self._cut_spectrogram(bg2, f1, f2)

    self.background=

def main():
    scriptname = 'ashortiso'
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='+', help='Name of the input file.\
')

    parser.add_argument('-t', '--itime', type=float, help='Injection time.')
    parser.add_argument('-td', '--tdeath', type=float, help='Time in which the isomer will have completly dissapeared after the injection time..')    
    parser.add_argument('-fc', '--fcen', type=float,  help='Frecuency center of the supposed isomer.')
    parser.add_argument('-fs', '--fspan', type=float, , help='Frecuency span around fcenter.')
    
    parser.add_argument('-v', '--verbose',
			help='Increase output verbosity', action='store_true')

    args = parser.parse_args()

    print(f'Running {scriptname} V{__version__}')
    if args.verbose: log.basicConfig(level=log.DEBUG)

    # here we go:                                                                                       
    log.info(f'File {args.filename} passed for processing the information of {args.refisotope}+{args.re\
fcharge}.')

if __name__ == '__main__':
    main()
