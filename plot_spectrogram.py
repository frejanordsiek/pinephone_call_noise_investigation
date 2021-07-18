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


__version__ = '0.2'


def plot_spectrogram(sound_file, csv_file, plot_file, nperseg,
                     focus_frequency=217.0, focus_window=4.0,
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
    focus_frequency : float, optional
        The frequency in Hz to center the last subplot on (the
        focus). Must be positive, finite, and non-NaN. The default is
        217 Hz.
    focus_window : float, optional
        The width in Hz of the frequency window to focus on in the last
        suplot. Must be positive, finite, and non-NaN. The default is 4
        Hz.
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
    if not isinstance(focus_frequency, float):
        raise TypeError('focus_frequency must be float.')
    if not numpy.isfinite(focus_frequency) or focus_frequency <= 0:
        raise ValueError('focus_frequency must be positive, non-NaN, '
                         'and finite.')
    if not isinstance(focus_window, float):
        raise TypeError('focus_window must be float.')
    if not numpy.isfinite(focus_window) or focus_window <= 0:
        raise ValueError('focus_window must be positive, non-NaN, '
                         'and finite.')
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

    # Make the plot. There will be four subplots, each above the other
    # that plot the spectrum with different fequency limits. The first
    # three gofrom 0 Hz up to an upper frequency limit
    # which starts from the maximum available and works its way down to
    # smaller and smaller fequencies. The fourth one zooms in on the
    # focus frequency.
    if plot_file is not None or showplot:
        fig, axs = plt.subplots(4, 1, figsize=(7.5, 10),
                                constrained_layout=True)
        for i, (lower_freq, upper_freq) in enumerate(
                ((None, None), (None, 5e3), (None, 500.0),
                 (focus_frequency - 0.5 * focus_window,
                  focus_frequency + 0.5 * focus_window))):
            ax = axs[i]
            if upper_freq is None:
                if lower_freq is None:
                    selector = slice(None)
                else:
                    selector = f >= lower_freq
            else:
                if lower_freq is None:
                    selector = f <= upper_freq
                else:
                    selector = numpy.logical_and(f >= lower_freq,
                                                 f <= upper_freq)
            ax.semilogy(f[selector], pwelch[selector], '-k')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Power Spectrum (A$^2$ / Hz)')
            if upper_freq is not None:
                if lower_freq is not None:
                    ax.set_xlim((lower_freq, upper_freq))
                else:
                    ax.set_xlim((0, upper_freq))
            elif lower_freq is not None:
                ax.set_xlim((lower_freq, ax.get_xlim()[1]))
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
    parser.add_argument('-n', type=int, default=2**18,
                        help='The number of samples per segment to use '
                        'when calculating the spectrogram. Must be '
                        'positve. Larger values increase the fequency '
                        'resolution but result in less segments to '
                        'average together (more noise). The default is '
                        '2**18 = 262144.')
    parser.add_argument('-f', '--focus-frequency', type=float,
                        default=217.0,
                        help='The center frequency in Hz to focus on '
                        'in the last subplot. The default is 217 Hz.')
    parser.add_argument('-w', '--focus-window', type=float,
                        default=4.0,
                        help='The full width of the fequency window '
                        'in Hz to focus on in the last subplot. The '
                        'default is 4 Hz.')
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
                     focus_frequency=args.focus_frequency,
                     focus_window=args.focus_window,
                     showplot=args.show_plot)


# Run main.
if __name__ == '__main__':
    sys.exit(main())
