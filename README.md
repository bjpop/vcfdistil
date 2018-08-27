[![travis](https://travis-ci.org/bjpop/vcfdistil.svg?branch=master)](https://travis-ci.org/bjpop/vcfdistil)

# Overview 

This program transforms and filters VCF files 

In the examples below, `$` indicates the command line prompt.

# Licence

This program is released as open source software under the terms of [BSD-3-Clause License](https://raw.githubusercontent.com/bjpop/vcfdistil/master/LICENSE).

# Installing

Clone this repository: 
```
$ git clone https://github.com/bjpop/vcfdistil
```

Move into the repository directory:
```
$ cd vcfdistil
```

Vcfdistil can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv vcfdistil_dev
$ source vcfdistil_dev/bin/activate
$ pip install -U /path/to/vcfdistil
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/vcfdistil
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/vcfdistil
```


# General behaviour

## Help message

Vcfdistil can display usage information on the command line via the `-h` or `--help` argument:

```
$ vcfdistil -h
```

## Logging

If the ``--log FILE`` command line argument is specified, vcfdistil will output a log file containing information about program progress. The log file includes the command line used to execute the program, and a note indicating which files have been processes so far. Events in the log file are annotated with their date and time of occurrence. 

```
$ vcfdistil --log bt.log file1.vcf file2.vcf
# normal vcfdistil output appears here
# contents of log file displayed below
```
```
$ cat bt.log
12/04/2016 19:14:47 program started
12/04/2016 19:14:47 command line: /usr/local/bin/vcfdistil --log bt.log file1.fasta file2.fasta
```

# Error handling

XXX FIXME

# Testing

## Unit tests

```
$ cd vcfdistil/python/vcfdistil
$ python -m unittest -v vcfdistil_test
```

## Test suite


# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[vcfdistil issue tracker](https://github.com/bjpop/vcfdistil/issues)
