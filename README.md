# ModifiedEMRIWaveforms
Repository for gravitational waves from asymmetric binaries with scalar fields, building on [FastEMRIWaveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms?tab=readme-ov-file#citation).

If you make use of this repository, please see the [citation section](#citation) below, together with the citation section in [FastEMRIWaveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms?tab=readme-ov-file#citation).

## Set up 

To make use of this repository and produce waveforms from asymmetric binaries with scalar fields, you need to install [FastEMRIWaveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms?tab=readme-ov-file#citation) first. A detailed installation guide can be found in the [official documentation](https://fastemriwaveforms.readthedocs.io/en/stable/). 


After installing FEW, clone this repository:

```
git clone https://github.com/susannabarsanti/ModifiedEMRIWaveforms.git
cd ModifiedEMRIWaveforms
pip install .
```
Alternatively, you can download the repository as a ZIP archive from GitHub.

## Repository Structure
The files in ModifiedEMRIWaveforms include: 
- mew/scalar_flux.py : python script containing the trajectory class 
- mew/data/ : data folder containing the fluxes for trajectory production
- notebooks/example.ipynb : interactive notebook demonstrating how to run the code

## Usage 
Once the repository is downloaded, you can produce your modified waveform in a python script by importing the class as
```
from mew import KerrCircEqFluxScalar
```
You can now use the trajectory class in FEW, for instance: 
```
from few.waveform import GenerateEMRIWaveform
inspiral_kwargs = {
    'flux_output_convention':'pex',
    'func':KerrCircEqFluxScalar
    }

Kerr_waveform = GenerateEMRIWaveform(
        "FastKerrEccentricEquatorialFlux",
        sum_kwargs=dict(pad_output=True),
        inspiral_kwargs = inspiral_kwargs,
        use_gpu=use_gpu,
        return_list=False,
    )
```

For further details, see notebooks/example.ipynb. 

## Software release

This repository is archived on Zenodo and each software release is assigned a DOI.

## Citation

If you use this code in your research, please cite the software release.
You can find the citation information in the `CITATION.cff` file in this repository.

## Acknowledgements

This code builds upon the
[FastEMRIWaveforms](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms)
framework developed by the Black Hole Perturbation Toolkit.

