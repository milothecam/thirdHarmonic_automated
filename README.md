# thirdHarmonic_automated

## What it does
These files are used to calculate the power at the input and output of a TKIP while measuring the generation of harmonics in it. The calculations are still to be verified and different approaches are being tested at the moment. To be updated.

## How to run
The script harmonic_measurement.py performs measurements using a signal generator and a VNA under Test Receiver mode (reading of signals with no output power from other ports of the VNA).
It will save csv files around the measured tones, to be used for further analysis.

The script harmonic_analysis_corrections_threeS2P.py is the current version of the analysis script. It plots the maximum power of all the saved csv files grouped by frequency from the signal generator tone and respective harmonic.

![image](C:\Users\cespinoz\Documents\Desktop\scriptsForWork\thirdHarmonic_automated\schematicNow.png)


## Requirements
- Python 3.x
- numpy
- etc.

## Example
Explain expected behavior or output.

## Notes
Anything important or weird.
