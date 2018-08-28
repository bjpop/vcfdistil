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
    parser.add_argument('--filter',
                        metavar='FILTER',
                        type=str,
                        help='file name of Python filter code')
    parser.add_argument('vcf_file',
                        nargs='?',
                        metavar='VCF_FILE',
                        type=str,
                        help='Input VCF file')
    return parser.parse_args()


def strip_quotes(text):
    if text.startswith('"') or text.startswith("'"):
        text = text[1:]
    if text.endswith('"') or text.endswith("'"):
        text = text[:-1]
    return text


def parse_metadata_record(record):
    if record.startswith('<') and record.endswith('>'):
        metadata_record_dict = {}
        no_brackets = record[1:-1]
        fields = no_brackets.split(',')
        for field in fields:
            key_value = field.split('=')
            if len(key_value) == 2:
                key, value = key_value[:2]
                unquoted_value = strip_quotes(value)
                metadata_record_dict[key] = unquoted_value 
        return metadata_record_dict 
    else:
        return None


metadata_regex = re.compile(r"\#\#(?P<name>[^=]*)=(?P<record>.*)$")


def parse_metadata(metadata_dict, row):
    match = re.match(metadata_regex, row)
    if match is not None:
        name = match.group('name')
        record = match.group('record')
        metadata_record = parse_metadata_record(record)
        if metadata_record is not None:
            if name not in metadata_dict:
                metadata_dict[name] = {}
            if 'ID' in metadata_record:
                this_id = metadata_record['ID']
                metadata_dict[name][this_id] = metadata_record
        else:
            metadata_dict[name] = record


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


class Record(object):
    '''
    info should be a dictionary
    genotypes, if present, should be a list
    '''
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format=None, genotypes=None):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.genotypes = genotypes 


    def __str__(self):
        fields = [self.chrom, self.pos, self.id, self.ref, self.alt, self.qual]
        if self.filter is not None:
            fields.append(self.filter)
        if self.genotypes is not None:
            fields.extend(self.genotypes)
        return '\t'.join(fields)

# INFO fields are encoded as a semicolon-separated series of short
# keys with optional values in the format: <key>=<data>[,data].

def parse_info(text):
    info_dict = {}
    fields = text.split(';')
    for field in fields:
        key_value = field.split('=')
        if len(key_value) == 2:
            key, value = key_value[:2]
            if ',' in value:
                value = value.split(',')
            info_dict[key] = value
    return info_dict



NUM_MANDATORY_RECORD_FIELDS = 8


def process_variants(file):
    for row in file:
        row = row.strip()
        fields = row.split('\t')
        if len(fields) >= NUM_MANDATORY_RECORD_FIELDS:
            chrom, pos, id, ref, alt, qual, filter, info = fields[:NUM_MANDATORY_RECORD_FIELDS]
            format = None
            genotypes = None
            if len(fields) >= NUM_MANDATORY_RECORD_FIELDS + 1:
                format = fields[NUM_MANDATORY_RECORD_FIELDS]
                genotypes = fields[NUM_MANDATORY_RECORD_FIELDS + 1:]
            info_dict = parse_info(info)
            yield Record(chrom, pos, id, ref, alt, qual, filter, info_dict, format, genotypes)


def process_vcf_file(file, begin, end, record_filter):
    metadata_iter, header_iter, variants_iter = tee(file, 3)
    metadata, num_metadata_rows = read_metadata(metadata_iter)
    sample_ids = read_header(islice(header_iter, num_metadata_rows, None, None))
    for item in begin(metadata, sample_ids):
        print(item)
    for record in process_variants(islice(variants_iter, num_metadata_rows + 1, None, None)):
        for filtered_record in record_filter(record, metadata, sample_ids):
            print(filtered_record)
    for item in end(metadata, sample_ids):
        print(item)


def process_file(options, begin, end, record_filter):
    '''
    '''
    if options.vcf_file is not None:
        logging.info("Processing VCF file from %s", options.vcf_file)
        try:
            vcf_file = open(options.vcf_file)
        except IOError as exception:
            exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
        else:
            with vcf_file:
                process_vcf_file(vcf_file, begin, end, record_filter)
    else:
        process_vcf_file(sys.stdin, begin, end, record_filter)


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


def identity_filter(x, _metadata, _sample_ids):
    yield x


def null_begin(_metadata, _sample_ids):
    return iter(())


def null_end(_metadata, _sample_ids):
    return iter(())


def get_custom_functions(options):
    # default filter is the identity function,
    # will include all rows in the output
    # WARNING: this code execs arbitrary Python code. Do not use this in
    # an untrusted environment, such as a web application!
    custom_filter = identity_filter 
    custom_begin = null_begin
    custom_end = null_end
    if options.filter:
        try:
            with open(options.filter) as custom_filter_file:
                this_locals = {}
                contents = custom_filter_file.read()
                exec(contents, globals(), this_locals)
                if 'filter' in this_locals:
                    custom_filter = this_locals['filter']
                if 'begin' in this_locals:
                    custom_begin = this_locals['begin']
                if 'end' in this_locals:
                    custom_end = this_locals['end']
        except Exception as e:
            print("Error when loading custom filter: {}".format(e))
            exit(1)
    return custom_begin, custom_end, custom_filter


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    begin, end, record_filter = get_custom_functions(options)
    init_logging(options.log)
    process_file(options, begin, end, record_filter)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
