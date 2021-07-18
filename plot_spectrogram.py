#!/usr/bin/env python3

# Copyright 2021 Freja Nordsiek
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import sys


import numpy
import scipy.io.wavfile
import scipy.signal
import matplotlib.pylab as plt


__version__ = '0.1'


def plot_spectrogram(sound_file, csv_file, plot_file, nperseg,
                     showplot=False):
    """ Plot the spectogram recorded from the phone's microphone.

    Parameters
    ----------
    sound_file : str
        The path to the sound file. Must be a WAV format.
    csv_file : str or None
        The path to write the spectrogram to in CSV form, or ``None`` to
        not write any. Can be optionally gzip compressed if it ends in a
        ``'.gz'`` extension.
    plot_file : str or None
        The path to write the plot of the spectrogram to, or ``None`` to
        not generate any plot. Must be an image format supported by
        ``matplotlib``, which is determined from the extension.
    nperseg : int
        The number of samples per segment to use when calculating the
        spectrogram. Must be positive.
    showplot : bool, optional
        Whether a plot of the spectrogram should be shown in a GUI
        window or not (default) at the end of processing.

    Raises
    ------
    TypeError
        If an argument has the wrong type.
    ValueError
        If an argument has an invalid value.
    IOError
        If a file cannot be read or writen.

    """
    # Check the arguments.
    if not isinstance(sound_file, str):
        raise TypeError('sound_file must be str.')
    if csv_file is not None and not isinstance(csv_file, str):
        raise TypeError('csv_file must be str or None.')
    if plot_file is not None and not isinstance(plot_file, str):
        raise TypeError('plot_file must be str or None.')
    if not isinstance(nperseg, int):
        raise TypeError('nperseg must be int.')
    if nperseg < 1:
        raise ValueError('nperseg must be positive.')
    if not isinstance(showplot, bool):
        raise TypeError('showplot must be bool')
    # Read the data and convert to float.
    sample_frequency, raw_data = scipy.io.wavfile.read(sound_file)
    data = raw_data.astype('float32')

    # Calculate the spectrogram for each channel and then average the
    # results across channels.
    f, pwelch_by_channel = scipy.signal.welch(data, sample_frequency,
                                              scaling='spectrum',
                                              nperseg=nperseg, axis=0)
    if pwelch_by_channel.ndim == 1:
        pwelch = pwelch_by_channel
    else:
        pwelch = pwelch_by_channel.mean(axis=1)

    # Write the spectrogram to the CSV file.
    if csv_file is not None:
        numpy.savetxt(csv_file, numpy.vstack((f, pwelch)).T,
                      fmt='%.9g', delimiter=', ', comments='# ',
                      header='Frequency (Hz), Power Spectrum (A^2/Hz)')

    # Make the plot. There will be three subplots, each above the other
    # that plot the spectrum from 0 Hz up to an upper frequency limit
    # which is different for each one. The upper frequency limits start
    # from the maximum available and work their way down to smaller and
    # smaller fequencies.
    if plot_file is not None or showplot:
        fig, axs = plt.subplots(3, 1, figsize=(7.5, 10),
                                constrained_layout=True)
        for i, upper_freq in enumerate((None, 5e3, 500.0)):
            ax = axs[i]
            if upper_freq is None:
                selector = slice(None)
            else:
                selector = f <= upper_freq
            ax.semilogy(f[selector], pwelch[selector], '-k')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Power Spectrum (A$^2$ / Hz)')
            if upper_freq is not None:
                ax.set_xlim((0, upper_freq))
            ax.grid('on')
            ax.tick_params(top=True, right=True)
        if plot_file is not None:
            fig.savefig(plot_file, dpi=150, pad_inches=0.0)
        if showplot:
            plt.show()


def _get_parser():
    """ Get the command line argument parser.

    Returns
    -------
    parser : argparse.ArgumentParser
        The argument parser.

    See Also
    --------
    argparse.ArgumentParser

    """
    parser = argparse.ArgumentParser(
        description='Plot the spectrogram of the audio recorded from '
        "the phone's microphone during a test call in order to see "
        "the noise pickup from the modem's transmitter.")

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + str(__version__))
    parser.add_argument('-c', '--csv', type=str,
                        default=None,
                        help='The output file path for the generated '
                        'spectrogram in numerical form. It will be a '
                        'CSV file, optionally gzip compressed based '
                        'on its extension. By default, no file is'
                        'written.')
    parser.add_argument('-p', '--plot', type=str,
                        default=None,
                        help='The output file path for the generated '
                        'plot of the spectrogram. It must be an image '
                        'format supported by matplotlib (based on the '
                        'extension). By default, no plot is generated.')
    parser.add_argument('-n', type=int, default=2**14,
                        help='The number of samples per segment to use '
                        'when calculating the spectrogram. Must be '
                        'positve. Larger values increase the fequency '
                        'resolution but result in less segments to '
                        'average together (more noise). The default is '
                        '2**14 = 16384.')
    parser.add_argument('-s', '--show-plot', action='store_true',
                        help='Show a plot of the spectrogram in a GUI '
                        'window when done processing.')
    parser.add_argument('soundfile', metavar='SOUNDFILE', type=str,
                        help='The sound file to read. Must be a WAV '
                        'file.')
    return parser


def main(argv=None):
    """ Run this program with it connected to the command line.

    Parameters
    ----------
    argv : iterable of str, optional
        The input arguments to use. The default is ``sys.argv``.

    See Also
    --------
    sys.argv

    """
    # Grab the input arguments if none were given.
    if argv is None:
        argv = sys.argv

    # Parse the input arguments.
    parser = _get_parser()
    args = parser.parse_args(argv[1:])

    # Run.
    plot_spectrogram(args.soundfile, args.csv, args.plot, args.n,
                     showplot=args.show_plot)


# Run main.
if __name__ == '__main__':
    sys.exit(main())
