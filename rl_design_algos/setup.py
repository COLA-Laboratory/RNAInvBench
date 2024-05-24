#create a setup.py file to install the package
from setuptools import setup, find_packages


setup(
    name="learna_tools",
    version="0.1.0",
    description="RNA Design with automated reinforcement learning.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        "Distance",
        "GraKeL==0.1.8",
        "numpy",
        "pandas",
        "pyaml",
        "pynisher",
        "pytest",
        "pylint",
        "psutil",
        "configspace",
        "scipy",
        "Sphinx",
        "tensorflow==1.4.0",
        "pynisher",
        "scikits.bootstrap",
        "hpbandster",
        "tqdm",
        "black",
        "Distance",
        "tensorforce==0.3.3",
        "dataclasses",
        "munch",
        "pyyaml",
        "scikit-learn",
        "LinearFold",
        "tabulate",
        "matplotlib",
        "biopython",
        "logomaker",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        # Add other relevant classifiers.
    ],
    scripts=[
        'bin/learna',
        'bin/meta-learna',
        'bin/meta-learna-adapt',
        'bin/liblearna',
    ],
# Optional
    author='Frederic Runge',
    author_email='runget@cs.uni-freiburg.de',
    license='MIT',
    keywords='RNA Design',
    url='https://github.com/Rungetf/learna_tools.git'
)