from . import utils
import os
import multiprocessing

"""
guess-encoding.py

The original guess-encoding.py script can be found in
<https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py>.

It was originally a software licensed under the MIT License reproduced bellow:

MIT License

Copyright (c) 2009-2011 Brent Pedersen, Haibao Tang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


"""
    awk 'NR % 4 == 0' your.fastq | python %prog [options]

guess the encoding of a stream of qual lines.
"""

encoding = {'Sanger': [33, (33, 73)], 'Solexa': [64, (59, 104)], 'Illumina-1.3': [64, (64, 104)], 'Illumina-1.5': [64, (66, 105)], 'Illumina-1.8': [33, (33, 74)]}


def get_qual_range(qual_str):
    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)


def get_encodings_in_range(rmin, rmax):
    valid_encodings = []
    for encoding_type, [phred, (emin, emax)] in encoding.items():
        if rmin >= emin and rmax <= emax:
            if encoding_type == 'Illumina-1.8':
                if rmax == emax:
                    valid_encodings.append(['Illumina-1.8', 33])
                else:
                    valid_encodings.append(['Sanger', 33])
            else:
                valid_encodings.append([encoding_type, phred])
    return valid_encodings if len(valid_encodings) > 0 else None


@utils.trace_unhandled_exceptions
def guess_encoding(fastq, number_reads_access_None_all, outdir):
    gmin, gmax = 99, 0
    valid_encodings = None
    reads_length = []
    with open(fastq, 'rtU') as reader:
        for i, line in enumerate(reader):
            if number_reads_access_None_all is None or (i + 1) / 4 <= number_reads_access_None_all:
                if (i + 1) % 4 == 0:
                    if len(line) > 0:
                        reads_length.append(len(line.splitlines()[0]))
                        lmin, lmax = get_qual_range(line.splitlines()[0])
                        if lmin < gmin or lmax > gmax:
                            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                            valid_encodings = get_encodings_in_range(gmin, gmax)

    utils.saveVariableToPickle([fastq, valid_encodings,
                                min(reads_length) if len(reads_length) > 0 else None,
                                max(reads_length) if len(reads_length) > 0 else None,
                                len(reads_length), sum(reads_length)],
                               outdir, 'encoding' + '.' + os.path.splitext(os.path.basename(fastq))[0])


def gather_data_together(data_directory):
    data = {}

    files = [f for f in os.listdir(data_directory) if not f.startswith('.') and os.path.isfile(os.path.join(data_directory, f))]
    for file_found in files:
        if file_found.startswith('encoding.') and file_found.endswith('.pkl'):
            file_path = os.path.join(data_directory, file_found)

            fastq, valid_encodings, min_reads_length, max_reads_length, num_reads, num_bp = \
                utils.extractVariableFromPickle(file_path)
            data[fastq] = {'valid_encodings': valid_encodings,
                           'min_reads_length': min_reads_length,
                           'max_reads_length': max_reads_length,
                           'num_reads': num_reads,
                           'num_bp': num_bp}

            os.remove(file_path)

    return data


def get_final_encoding(encoding_data):
    possible_encoding = {}
    for fastq in encoding_data:
        if encoding_data[fastq] is not None and encoding_data[fastq]['valid_encodings'] is not None:
            for i in range(0, len(encoding_data[fastq]['valid_encodings'])):
                if encoding_data[fastq]['valid_encodings'][i][0] not in possible_encoding:
                    possible_encoding[encoding_data[fastq]['valid_encodings'][i][0]] = 0
                possible_encoding[encoding_data[fastq]['valid_encodings'][i][0]] += 1

    final_encoding = []
    if len(possible_encoding) > 0:
        if list(possible_encoding.values()).count(max(possible_encoding.values())) > 0:
            for encoding_type in possible_encoding:
                if possible_encoding[encoding_type] == max(possible_encoding.values()):
                    final_encoding = [encoding_type, encoding[encoding_type][0]]
        else:
            final_encoding = None
    else:
        final_encoding = None

    return final_encoding


def determine_min_max_reads_length(encoding_data):
    """
    Returns the minimum and maximum reads length found for all fastq and for each fastq

    Parameters
    ----------
    encoding_data : dict
        Dictionary with encondig data, and reads length for each fastq. Something like
        data[fastq] = {'valid_encodings': valid_encodings,
                       'min_reads_length': min_reads_length,
                       'max_reads_length': max_reads_length,
                       'num_reads': num_reads,
                       'num_bp': num_bp}

    Returns
    -------
    min_reads_length_found : int
        Minimum reads length found in fastq file
    max_reads_length_found : int
        Maximum reads length found in fastq file
    min_reads_length_each_fastq : list
        Minimum reads length found for each fastq file
    max_reads_length_each_fastq : list
        Maximum reads length found for each fastq file
    """
    min_length_each_fastq = [encoding_data[fastq]['min_reads_length'] for fastq in encoding_data if
                             encoding_data[fastq]['min_reads_length'] is not None]

    max_length_each_fastq = [encoding_data[fastq]['max_reads_length'] for fastq in encoding_data if
                             encoding_data[fastq]['max_reads_length'] is not None]

    return min(min_length_each_fastq) if len(min_length_each_fastq) > 0 else None, \
           max(max_length_each_fastq) if len(max_length_each_fastq) > 0 else None, \
           min_length_each_fastq, \
           max_length_each_fastq


def get_num_reads_bp(encoding_data):
    """
    Returns the total number of reads and bp sequenced

    Parameters
    ----------
    encoding_data : dict
        Dictionary with encondig data, and reads length for each fastq. Something like
        data[fastq] = {'valid_encodings': valid_encodings,
                       'min_reads_length': min_reads_length,
                       'max_reads_length': max_reads_length,
                       'num_reads': num_reads,
                       'num_bp': num_bp}

    Returns
    -------
    num_reads : int
        Total number of reads sequenced
    num_bp : int
        Total number of bp sequenced
    """
    num_reads = [encoding_data[fastq]['num_reads'] for fastq in encoding_data if
                 encoding_data[fastq]['num_reads'] is not None]

    num_bp = [encoding_data[fastq]['num_bp'] for fastq in encoding_data if
              encoding_data[fastq]['num_bp'] is not None]

    return sum(num_reads) if len(num_reads) > 0 else None, sum(num_bp) if len(num_bp) > 0 else None


def fastq_files_enconding(fastq_files_list, number_reads_access_None_all, outdir, threads):
    pool = multiprocessing.Pool(processes=threads)
    for fastq in fastq_files_list:
        pool.apply_async(guess_encoding, args=(fastq, number_reads_access_None_all, outdir,))
    pool.close()
    pool.join()

    encoding_data = gather_data_together(outdir)

    final_encoding = get_final_encoding(encoding_data)

    min_reads_length, max_reads_length, _, _ = determine_min_max_reads_length(encoding_data)

    num_reads, num_bp = get_num_reads_bp(encoding_data)

    return final_encoding, min_reads_length, max_reads_length, num_reads, num_bp
