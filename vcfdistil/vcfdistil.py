'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 27 Aug 2018 
License     : BSD-3-Clause 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

'''

from argparse import ArgumentParser
import sys
import logging
from itertools import takewhile, tee, islice
import re
import pkg_resources


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_VCF_FILE_ERROR = 3
DEFAULT_VERBOSE = False
PROGRAM_NAME = "vcfdistil"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more VCF files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('vcf_files',
                        nargs='*',
                        metavar='VCF_FILE',
                        type=str,
                        help='Input VCF files')
    return parser.parse_args()


metadata_regex = re.compile(r"\#\#(?P<key>[^=]*)=(?P<value>.*)$")


def parse_metadata(metadata_dict, row):
    match = re.match(metadata_regex, row)
    if match is not None:
        key = match.group('key')
        value = match.group('value')
        metadata_dict[key] = value


def read_metadata(file):
    metadata_dict = {}
    num_metadata_rows = 0
    for metadata_row in takewhile(lambda x: x.startswith('##'), file):
        parse_metadata(metadata_dict, metadata_row)
        num_metadata_rows += 1
    return metadata_dict, num_metadata_rows


mandatory_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


def parse_header(row):
    '''returns the sample ids, if they are present'''
    fields = row.split('\t')
    num_mandatory = len(mandatory_header_fields)
    print(fields[:num_mandatory])
    if fields[:num_mandatory] == mandatory_header_fields:
        remaining_fields = fields[num_mandatory:]
        if len(remaining_fields) > 0 and remaining_fields[0] == "FORMAT":
            sample_ids = remaining_fields[1:]
            return sample_ids
    else:
        logging.warning("header missing mandatory fields: {}".format(row))
    return None


def read_header(file):
    for header_row in takewhile(lambda x: x.startswith('#CHROM'), file):
        return parse_header(header_row.strip())
    return None


def process_variants(file):
    for row in file:
        row = row.split()
        print(row)


def process_vcf_file(file):
    metadata_iter, header_iter, variants_iter = tee(file, 3)
    metadata, num_metadata_rows = read_metadata(metadata_iter)
    print((metadata, num_metadata_rows))
    sample_ids = read_header(islice(header_iter, num_metadata_rows, None, None))
    print(sample_ids)
    process_variants(islice(variants_iter, num_metadata_rows + 1, None, None))


def process_files(options):
    '''
    '''
    for vcf_filename in options.vcf_files:
        logging.info("Processing VCF file from %s", vcf_filename)
        try:
            vcf_file = open(vcf_filename)
        except IOError as exception:
            exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
        else:
            with vcf_file:
                process_vcf_file(vcf_file)


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
