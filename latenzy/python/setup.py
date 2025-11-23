from setuptools import setup, find_packages

setup(
    name="latenzy",
    version="0.1.0",
    description="Python tools for latenZy analysis",
    long_description="latenzy.py\n\n"
                     "Contains latenzy and latenzy2 to compute latencies for spiking responses\n"
                     "See Haak et al. 2025\n\n"
                     "2025, Robin Haak, Alexander Heimel",
    author="Robin Haak, Alexander Heimel",
    url="https://github.com/Herseninstituut/latenZy",
    license="GNU General Public License v3 (GPL-3.0)",
    packages=find_packages('.'),
    package_dir={'': '.'},
    install_requires=[],
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Applications",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    include_package_data=True,
    zip_safe=False,
)
