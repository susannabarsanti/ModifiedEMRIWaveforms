from setuptools import setup, find_packages

setup(
    name="ModifiedEMRIWaveform",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "h5py",
    ],
    include_package_data=True,
    package_data={
        "mew": ["data/*.dat"],
    },
)
