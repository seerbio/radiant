# -*- coding: utf-8 -*-
"""Test program using Python wrapper for timsdata.dll to read tsf (spectrum) data"""

import sys, tsfdata, sqlite3

if len(sys.argv) < 2:
    raise RuntimeError("need arguments: tsf_directory")

analysis_dir = sys.argv[1]

if sys.version_info.major == 2:
    analysis_dir = unicode(analysis_dir)

tsf = tsfdata.TsfData(analysis_dir)
conn = tsf.conn


# Get total spectrum count:
q = conn.execute("SELECT COUNT(*) FROM Frames")
row = q.fetchone()
N = row[0]
print("Analysis has {0} spectra.".format(N))

# What type of spectra are contained
has_line_spectra = conn.execute("SELECT Value FROM GlobalMetadata WHERE Key='HasLineSpectra'").fetchone()[0] == '1'
has_line_spectra_peak_width = conn.execute("SELECT Value FROM GlobalMetadata WHERE Key='HasLineSpectraPeakWidth'").fetchone()[0] == '1'
has_profile_spectra = conn.execute("SELECT Value FROM GlobalMetadata WHERE Key='HasProfileSpectra'").fetchone()[0] == '1'

print(f"Analysis has {'' if has_line_spectra else 'no '}line spectra {'with peak width information ' if has_line_spectra_peak_width else ''}and {'' if has_profile_spectra else 'no '}profile spectra.")

spectrum_id = 1

# Get a line spectrum:
if has_line_spectra:
    if has_line_spectra_peak_width:
        indices, intensities, widths = tsf.readLineSpectrumWithWidth(spectrum_id)
        peak_mz = tsf.indexToMz(spectrum_id, indices)
        print("Line spectrum with width:")
        print(peak_mz)
        print(intensities)
        print(widths)
    else:
        indices, intensities = tsf.readLineSpectrum(spectrum_id)
        peak_mz = tsf.indexToMz(spectrum_id, indices)
        print("Line spectrum without width:")
        print(peak_mz)
        print(intensities)

# Get a profile spectrum
if has_profile_spectra:
    intensities = tsf.readProfileSpectrum(spectrum_id)
    indices = list(range(len(intensities)))
    mz = tsf.indexToMz(spectrum_id, indices)
    print("Profile spectrum:")
    print(mz)
    print(intensities)
