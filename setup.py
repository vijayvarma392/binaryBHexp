'''
Standard setup.py to upload the code on pypi.
Do:
    python setup.py sdist bdist_wheel
    twine upload dist/*
'''
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# Extract code version from __init__.py
def get_version():
    with open('binaryBHexp') as f:
        for line in f.readlines():
            if "__version__" in line:
                return line.split('"')[1]

setuptools.setup(
    name="binaryBHexp",
    version=get_version(),
    author="Vijay Varma",
    author_email="vvarma@caltech.edu",
    description="The binary Black Hole explorer (visualizations).",
    keywords='black-holes gravitational-waves visualization',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vijayvarma392/binaryBHexp",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'palettable',
        'surfinBH',
        'NRSur7dq2>=1.0.5',
    ],
    scripts=['binaryBHexp'],
    classifiers=[
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
    ],
)
