#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''The transforms and filters VCF files.'''


setup(
    name='vcfdistil',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['vcfdistil'],
    package_dir={'vcfdistil': 'vcfdistil'},
    entry_points={
        'console_scripts': ['vcfdistil = vcfdistil.vcfdistil:main']
    },
    url='https://github.com/bjpop/vcfdistil',
    license='LICENSE',
    description=('Transform and filter VCF files'),
    long_description=(LONG_DESCRIPTION),
    install_requires=[],
)
