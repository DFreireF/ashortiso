import argparse
import logging as log
#from ashortiso.version import __version__
from iqtools import *
from lmfit import *
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


class IsomerIdentification():
    def __init__(self, filename,  initial_time, time_length, fcen, fspan, lframes):
        self.fcen = fcen
        self.fspan = fspan
        self.itime = initial_time
        self.tdeath = time_length
        self.filename = filename
        self.lframes = lframes

    def _create_spectrogram(self, skip, method=None): #i like this one
        iq = get_iq_object(self.filename)
        iq.read_samples(1)
        nframes = int(self.tdeath*iq.fs/self.lframes)
        sframes = int(skip*iq.fs/self.lframes)
        iq.read(nframes=nframes, lframes=self.lframes, sframes=sframes)
        iq.method = 'mtm' #'fft', 'mtm', 'welch'
        if method: iq.method = method
        return iq.get_spectrogram(nframes, self.lframes)
    
    def _get_averaged_window(self, skip, fcen=0, fspan=1.5e3, plot=False, savefig=False):
        xx, yy, zz = self._create_spectrogram(skip)
        nxx, nyy, nzz = get_cut_spectrogram(xx, yy, zz,xcen=0, xspan=fspan)
        axx, ayy, azz = get_averaged_spectrogram(nxx, nyy, nzz, len(nxx[:,0]))
        if plot:
            fig, axs = plt.subplots(3,1)
            fig.suptitle('Spectrograms + Averaged spectrum')
            axs[0].plot_spectrogram(xx, yy, zz)
            axs[1].plot_spectrogram(nxx, nyy, nzz)
            axs[2].plot_spectrum(axx[0,:], azz[0,:], dbm=True)
            if savefig: plt.savefig(datetime.now().strftime('%Y.%m.%d_%H.%M.%S')+'.plot.spectrums.pdf')
        return axx[0,:], azz[0,:]

    def get_isomer_window(self, skip, fspan, fcen=0, fiso=-2e3): #fiso respect mother
        x, y = self._get_averaged_window(skip, fcen, fspan*10, plot=True)
        nxcen, area_gs = IsomerIdentification.fit_gaussian(x,y, amp=6e-07, cen=2e3, wid=2e2)
        xi, yi = self._get_averaged_window(skip, nxcen-fiso, fspan, plot=True)
        return xi, yi, area_gs

    def method_1(self, fspan, fcen=0, tspan=3):
        xi, yi, area_gsi = self.get_isomer_window(self.itime, fspan, fcen)
        energy_isomer_i = IsomerIdentification.energy_in_window_discrete(xi, yi)
        xf, yf, area_gsf = self.get_isomer_window(self.itime+tspan, fspan, fcen)
        energy_isomer_f = IsomerIdentification.energy_in_window_discrete(xf, yf)
        isomer=IsomerIdentification.isomer_or_not(energy_isomer_i, energy_isomer_f, factor=energy_isomer_i/energy_isomer_f)
        return isomer
    
    @staticmethod
    def evolution(x, y, z, plot=False, show=True): #x=frec, y=time, z=area (energy)
        delta_x, delta_y, delta_z= [np.array([]) for i in range(0,3)]
        for i in range(len(x)-1):
            delta_x=np.append(delta_x, x[i+1]-x[i])
            delta_y=np.append(delta_y, (y[i+1]-y[i]) + y[i])
            delta_e=np.append(delta_e, z[i+1]-z[i])
        if plot:
            fig, axs = plt.subplots(2,1, sharex=True)
            fig.suptitle('Variation of energy and frecuency of a peak with time')
            axs[0].plot(delta_y, delta_z)
            axs[1].plot(delta_y, delta_x)
            if show: plt.show()
            else: plt.savefig(datetime.now().strftime('%Y.%m.%d_%H.%M.%S')+'.plot.deltas-time.pdf')

    @staticmethod
    def isomer_or_not(isomer_e, back_e, factor=1):
        delta_e = np.abs(isomer_e-back_e)
        if delta_e > factor*back_e: return True
        else: return False

    @staticmethod
    def gaussian(x, amp, cen, wid):
        """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return (amp/(np.sqrt(2*np.pi)*wid))*np.exp(-(x-cen)**2/(2*wid**2))

    @staticmethod
    def fit_gaussian(x, y, amp, cen, wid, plot=False, fit_report=False, savefig=False):
        gmod = Model(IsomerIdentification.gaussian)
        result = gmod.fit(y, x = x, amp = 6e-07, cen = 0, wid = 2e2)
        new_xcen = result.params['cen'].value
        area = result.params['amp'].value
        if plot:
            plt.plot(x, y, 'o', label='data')
            plt.plot(x, result.init_fit, '--', label='initial fit')
            plt.plot(x, result.best_fit, '-', label='best fit')
            plt.show()
        if fit_report: print(result.fit_report())
        if savefig: plt.savefig(datetime.now().strftime('%Y.%m.%d_%H.%M.%S')+'.plot.gaussian.pdf')
        return new_xcen, area

    @staticmethod    
    def energy_in_window_discrete(x, y):
        delta_x = x[1]-x[0]
        return np.sum(y)*delta_x

    @staticmethod    
    def energy_in_window_cont(x, y):
        from scipy import integrate
        return integrate.simpson(y, x)


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

def files_with_isomers_or_not(file, itime, timetostudy, fcen, fspan, binning, fwithiso, fwithoutiso, wf):
    iso=IsomerIdentification(file, itime, timetostudy, fcen, fspan, binning)
    isomer=iso.method_1(fspan)
    if isomer:
        np.savetxt(wf,['##ISO##'+file+'##ISO##'], newline='\n', fmt='%s')
        fwithiso = fwithiso+1
    else:
        np.savetxt(wf, [file], newline='\n', fmt='%s')
        fwithoutiso = fwithoutiso+1
    return fwithiso, fwithoutiso

def print_output(fwith,fwithout):
    total_files_analysed=fwith+fwithout
    print(f'It has been analysed {total_files_analysed} with isomers present in {fwith} files and {fwithout} files without evidence of isomers.')

if __name__ == '__main__':
    main()
