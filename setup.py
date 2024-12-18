from setuptools import setup #, find_packages
#from codecs import open
#from os import path

#here = path.abspath(path.dirname(__file__))

setup(
    name='clumping',

    version='0.1',

    author='Phil Bradley and Stefan Schattgen',
    author_email= 'pbradley@fredhutch.org and Stefan.Schattgen@stjude.org',

    license='MIT',

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
    ],

    install_requires=["pandas", "numpy", "scipy","scikit-learn"],

    packages=['clumping','clumping.tcrdist'],
    #packages=find_packages(),

    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)
