Overview
========

During phone calls with poor signal with my PinePhone, the person on the other
end noticed a loud noise coming through that made it very hard to hear me. On my
end, I heard a sound coming from the inside of the phone that sounded like my
past phones when their transmit power was set to high. We managed to get it
recorded on a recorder and thus the investigation into this sound. My progress
can be read on the Pine64 Forum
(https://forum.pine64.org/showthread.php?tid=13922).

The tools here to quantitatively investigate the sound. On Mobian (and likely
other distros), the sound system is configured during a call in a way that it is
still possible to record the sound from the builtin microphone while the call is
going with the standard Sound Recorder app. The sound can be recorded during the
call in a quite place, exported to a FLAC file by the app, transfered to a
computer, and converted to a WAV file. The tools here are for analyzing the
results.

Everything in this repository is licensed under the MIT license
(https://github.com/frejanordsiek/blob/main/COPYING.txt).


Making Spectrograms of The Noise
================================

:file:`plot_spectrogram.py` is a script for taking a recorded WAV file,
computing the spectrogram, and exporting the spectrogram to a CSV file
(optionally gzip compressed) and/or as a plot (any image format supported by
`matplotlib <https://matplotlib.org>`_).

The ``-c CSVFILE`` option specifies the name of the CSV file to write the data
to.

The ``-p PLOTFILE`` option specifies the name of the file to write the plot
to. The extension determines the image format.

The ``-n NPERSEG`` option specifies the number of samples per segment to use
when computing the spectrogram. The default is ``2**18 = 262144`` (6 seconds at
44 kHz), which is a good value to be able to resolve small frequency peaks
while averaging enough segments together to get a reasonably smooth spectrogram
for a recording of about a minute.

The ``-f FOCUSFREQUENCY`` and ``-w FOCUSWINDOW`` options specify a frequency
window (center and width in Hz) to focus on for the last subplot of the
spectrogram plot in order to focus on a particular region in the spectrogram.
The default, 217 Hz with a width of 4 Hz, is centered on the first peak
of the fequency comb on one particular pinephone.

The ``-s`` option causes the spectrogram to be plotted in a GUI window when
done to allow for easier viewing than a plot file provides.

The one mandatory argument (the last one) is the path to the sound file to read
(WAV format).

I am keeping recorded spectrograms in the :file:`spectrograms` (CSV data) and
:file:`plots` directories. Feel free to add any with a Pull Request.
