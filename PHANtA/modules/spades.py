from . import utils
import os
import shutil
from functools import partial
import time
import re


# Run Spades
def spades(spades_folder, threads, fastq_files, notUseCareful, maxMemory, minCoverageAssembly, kmers,
           assembled_se_reads, spades_version=None, estimated_coverage=None, spades_not_use_isolate=False):
    contigs = os.path.join(spades_folder, 'contigs.fasta')

    command = ['spades.py', '', '--only-assembler', '--threads', str(threads), '--memory', str(maxMemory),
               '--cov-cutoff', str(minCoverageAssembly), '', '-1', fastq_files[0], '-2', fastq_files[1], '',
               '-o', spades_folder]

    if not notUseCareful:
        command[1] = '--careful'

    if not spades_not_use_isolate and \
            spades_version is not None and float('.'.join(spades_version.split('.')[:2])) >= 3.14 and \
            estimated_coverage is not None and estimated_coverage >= 100:
        command[1] = '--isolate'  # Isolate mode already implies --only-assembler, so this option has no effect.

    if len(kmers) > 0:
        kmers = ','.join(map(str, kmers))
        command[9] = str('-k ' + kmers)

    if assembled_se_reads is not None:
        command[14] = str('-s ' + assembled_se_reads)

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully, contigs


def define_kmers(kmers, maximumReadsLength):
    kmers_use = []
    if maximumReadsLength is not None:
        if kmers is None:
            if maximumReadsLength >= 175:
                kmers = [55, 77, 99, 113, 127]
            else:
                kmers = [21, 33, 55, 67, 77]
        for kmer in kmers:
            if kmer <= maximumReadsLength:
                kmers_use.append(kmer)
    return kmers_use


def define_minContigsLength(maximumReadsLength, minContigsLength):
    minimum_length = 200
    if minContigsLength is not None:
        minimum_length = minContigsLength
    else:
        if maximumReadsLength > minimum_length:
            minimum_length = maximumReadsLength

    return minimum_length


def define_memory(maxMemory, threads, available_memory_GB, spades_version=None):
    GB_per_thread = 2048 / 1024.0
    if spades_version is not None and float('.'.join(spades_version.split('.')[:2])) >= 3.14:
        GB_per_thread = 3072 / 1024.0

    minimum_required_memory_GB = GB_per_thread * threads
    if spades_version is not None and float('.'.join(spades_version.split('.')[:2])) >= 3.14:
        if minimum_required_memory_GB < 6:
            minimum_required_memory_GB = 6
    else:
        if minimum_required_memory_GB < 4:
            minimum_required_memory_GB = 4

    if available_memory_GB == 0:
        print('WARNING: it was not possible to determine the free available memory!')
        if maxMemory is not None:
            print('Setting SPAdes maximum memory to the one provided by user')
            time.sleep(10)
            return maxMemory
        else:
            print('Trying use the minimum memory required for SPAdes to run')
            available_memory_GB = minimum_required_memory_GB

    if maxMemory is None:
        if minimum_required_memory_GB > available_memory_GB:
            print('WARNNING: the minimum memory required to run SPAdes with {threads} threads ({min_mem} GB) are higher'
                  ' than the available memory ({avail_mem} GB)!'.format(threads=threads,
                                                                        min_mem=round(minimum_required_memory_GB, 1),
                                                                        avail_mem=round(available_memory_GB, 1)))
            print('Setting SPAdes maximum memory to {max_mem} GB'.format(max_mem=round((available_memory_GB - 0.5), 0)))
            return int(round((available_memory_GB - 0.5), 0))
        else:
            print('Setting SPAdes maximum memory to {min_mem} GB'.format(min_mem=round(minimum_required_memory_GB, 0)))
            return int(round(minimum_required_memory_GB, 0))
    else:
        if maxMemory < minimum_required_memory_GB < available_memory_GB:
            print('WARNNING: the minimum memory required to run SPAdes with {threads} threads ({min_mem} GB) are higher'
                  ' than the maximum memory set ({max_mem} GB)!'.format(threads=threads,
                                                                        min_mem=round(minimum_required_memory_GB, 1),
                                                                        max_mem=maxMemory))
            print('Consider to increase the SPAdes maximum memory'
                  ' to {min_mem} GB'.format(min_mem=round(minimum_required_memory_GB, 0)))
        elif maxMemory > available_memory_GB:
            print('WARNNING: the maximum memory set are higher than the available'
                  ' memory ({avail_mem} GB)!'.format(avail_mem=round(available_memory_GB, 1)))
            if minimum_required_memory_GB > available_memory_GB:
                print('WARNNING: the minimum memory required to run SPAdes with {threads} threads ({min_mem} GB) are'
                      ' higher than the available'
                      ' memory ({avail_mem} GB)!'.format(threads=threads, min_mem=round(minimum_required_memory_GB, 1),
                                                         avail_mem=round(available_memory_GB, 1)))
                print('Nevertheless, consider setting SPAdes maximum memory'
                      ' to {max_mem} GB'.format(max_mem=round((available_memory_GB - 0.5), 0)))
            else:
                print('Consider setting SPAdes maximum memory to {max_mem} GB (the minimum memory required to run'
                      ' SPAdes with {threads} threads)'.format(max_mem=round(minimum_required_memory_GB, 0),
                                                               threads=threads))
        time.sleep(10)
        return maxMemory


def get_SPAdes_sequence_information(spadesContigs):
    sequence_dict = {}

    with open(spadesContigs, 'rtU') as original_sequences:
        blank_line_found = False
        sequence_counter = 0
        for line in original_sequences:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not blank_line_found:
                    if line.startswith('>'):
                        sequence_counter += 1
                        sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': [], 'length': 0, 'AT': 0, 'GC': 0, 'N': 0, 'kmer_cov': float(line.split('_')[5]), 'discard': True}
                    else:
                        sequence_dict[sequence_counter]['sequence'].append(line)
                        sequence_dict[sequence_counter]['length'] += len(line)
                        line = line.upper()
                        sequence_dict[sequence_counter]['AT'] += len(re.findall('[AT]', line))
                        sequence_dict[sequence_counter]['GC'] += len(re.findall('[GC]', line))
                        sequence_dict[sequence_counter]['N'] += len(re.findall('[^ATGC]', line))
                else:
                    sequence_dict = None
            else:
                blank_line_found = True

    return sequence_dict


def determine_sequences_to_filter(sequence_dict, minContigsLength, minCoverageContigs, min_GC_content):
    for i in sequence_dict:
        sequence_dict[i]['discard'] = True

    spades_report_general = {'original': {'contigs': 0, 'bp': 0}, 'filtered': {'contigs': 0, 'bp': 0}}

    spades_report_general['original']['contigs'] = len(sequence_dict)
    spades_report_general['original']['bp'] = sum(sequence_dict[i]['length'] for i in sequence_dict)

    for i in sequence_dict:
        if sequence_dict[i]['length'] >= minContigsLength and \
                sequence_dict[i]['kmer_cov'] >= minCoverageContigs and \
                (sequence_dict[i]['GC'] / float(sequence_dict[i]['length']) >= min_GC_content and
                 sequence_dict[i]['GC'] / float(sequence_dict[i]['length']) <= 1 - min_GC_content):
            sequence_dict[i]['discard'] = False
            spades_report_general['filtered']['contigs'] += 1
            spades_report_general['filtered']['bp'] += sequence_dict[i]['length']

    return sequence_dict, spades_report_general


def write_filtered_sequences_and_stats(sequence_dict, spades_report_general, original_sequence_file, filtered_sequence_file, sampleName, write_only_report_original_True, saveExcludedContigs):
    if write_only_report_original_True is False:
        report_filtered = open(os.path.join(os.path.dirname(filtered_sequence_file), str('spades_report.filtered.' + os.path.splitext(os.path.basename(filtered_sequence_file))[0].split('.', 1)[1]) + '.tab'), 'wt')

    if saveExcludedContigs:
        path_excluded_contigs = os.path.splitext(filtered_sequence_file)[0] + '.excluded_contigs.fasta'
        excluded_contigs = open(path_excluded_contigs, 'wt')

    found_excluded_contigs = False
    with open(os.path.join(os.path.dirname(filtered_sequence_file), str('spades_report.original.' + os.path.splitext(os.path.basename(original_sequence_file))[0].split('.', 1)[1]) + '.tab'), 'wt') as report_original:
        with open(filtered_sequence_file, 'wt') as contigs_filtered:
            fields = ['header', 'length', 'AT', 'GC', 'N', 'kmer_cov']

            report_original.write('\n'.join(['#general', '>contigs', str(spades_report_general['original']['contigs']), '>bp', str(spades_report_general['original']['bp'])]) + '\n')
            report_original.write('#' + '\t'.join(fields) + '\n')
            if write_only_report_original_True is False:
                report_filtered.write('\n'.join(['#general', '>contigs', str(spades_report_general['filtered']['contigs']), '>bp', str(spades_report_general['filtered']['bp'])]) + '\n')
                report_filtered.write('#' + '\t'.join(fields) + '\n')

            for i in range(1, len(sequence_dict) + 1):
                report_original.write('\t'.join([str(sequence_dict[i][f]) for f in fields]) + '\n')
                if not sequence_dict[i]['discard']:
                    contigs_filtered.write('>' + sampleName + '_' + sequence_dict[i]['header'] + '\n' + '\n'.join(sequence_dict[i]['sequence']) + '\n')
                    if write_only_report_original_True is False:
                        report_filtered.write('\t'.join([str(sequence_dict[i][f]) for f in fields]) + '\n')
                else:
                    if saveExcludedContigs:
                        found_excluded_contigs = True
                        excluded_contigs.write('>' + sampleName + '_' + sequence_dict[i]['header'] + '\n' + '\n'.join(sequence_dict[i]['sequence']) + '\n')

    if saveExcludedContigs:
        excluded_contigs.flush()
        excluded_contigs.close()
        if not found_excluded_contigs:
            os.remove(path_excluded_contigs)


def qc_assembly(spades_report_general, estimatedGenomeSizeMb, maxNumberContigs):
    warnings = {}
    warnings['sample'] = False
    minimumBP = False

    if spades_report_general['filtered']['bp'] < estimatedGenomeSizeMb * 1000000 * 0.8 or spades_report_general['filtered']['bp'] > estimatedGenomeSizeMb * 1000000 * 1.5:
        warnings['sample'] = 'The number of assembled nucleotides (' + str(spades_report_general['filtered']['bp']) + ') are lower than 80% or higher than 150% of the provided estimated genome size'
    else:
        maxContigs = estimatedGenomeSizeMb * maxNumberContigs / 1.5
        if spades_report_general['filtered']['contigs'] > maxContigs:
            warnings['sample'] = 'The number of assembled contigs ({contigs}) exceeds {maxContigs}'.format(contigs=spades_report_general['filtered']['contigs'], maxContigs=maxContigs)
            print(warnings['sample'])

    if spades_report_general['filtered']['bp'] >= estimatedGenomeSizeMb * 1000000 * 0.8:
        minimumBP = True

    return warnings, minimumBP


def decide_filter_parameters(sequence_dict, minContigsLength, minCoverageContigs, estimatedGenomeSizeMb, maxNumberContigs):
    warnings = {}
    warnings['sample'] = False

    min_GC_content = 0.05

    filtered_sequences_sufix = 'length_GCcontent_kmerCov'

    print('Filtering for contigs with at least {min_len} nucleotides, a k-mer coverage of {min_cov} and a CG content'
          ' between {min_gc}% and {max_gc}%'.format(min_len=minContigsLength, min_cov=minCoverageContigs,
                                                    min_gc=min_GC_content * 100, max_gc=(1 - min_GC_content) * 100))
    sequence_dict, spades_report_general = determine_sequences_to_filter(sequence_dict, minContigsLength,
                                                                         minCoverageContigs, min_GC_content)

    warnings, minimumBP = qc_assembly(spades_report_general, estimatedGenomeSizeMb, maxNumberContigs)
    if warnings['sample'] is not False and not minimumBP:
        print('WARNING: ' + str(spades_report_general['filtered']['bp']) + ' assembled nucleotides in ' + str(spades_report_general['filtered']['contigs']) + ' contigs')

        filtered_sequences_sufix = 'length_GCcontent'

        print('Filtering for contigs with at least ' + str(minContigsLength) + ' nucleotides and a CG content between ' + str(min_GC_content * 100) + '% and ' + str((1 - min_GC_content) * 100) + '%')
        sequence_dict, spades_report_general = determine_sequences_to_filter(sequence_dict, minContigsLength, 0, min_GC_content)

        warnings, minimumBP = qc_assembly(spades_report_general, estimatedGenomeSizeMb, maxNumberContigs)
        if warnings['sample'] is not False and not minimumBP:
            print('WARNING: ' + str(spades_report_general['filtered']['bp']) + ' assembled nucleotides in ' + str(spades_report_general['filtered']['contigs']) + ' contigs')

    if warnings['sample'] is False:
        warnings = {}

    return warnings, sequence_dict, filtered_sequences_sufix, spades_report_general


spades_timer = partial(utils.timer, name='SPAdes')


# Run SPAdes procedure
@spades_timer
def run_spades(sample_name, outdir, threads, fastq_files, not_use_careful, max_memory, min_coverage_assembly,
               min_contigs_length, estimated_genome_size_mb, kmers, maximum_reads_length, default_kmers,
               min_coverage_contigs, assembled_se_reads, save_excluded_contigs, max_number_contigs,
               keep_scaffolds=False, spades_version=None, estimated_coverage=None, spades_not_use_isolate=False):
    pass_qc = True
    failing = {'sample': False}
    warnings = {}

    # Create SPAdes output directory
    spades_folder = os.path.join(outdir, 'spades', '')
    utils.removeDirectory(spades_folder)
    os.mkdir(spades_folder)

    # Determine k-mers to run
    if default_kmers:
        kmers = []
    else:
        kmers = define_kmers(kmers, maximum_reads_length)
        if len(kmers) == 0:
            print('SPAdes will use its default k-mers')
        else:
            print('SPAdes will use the following k-mers: ' + str(kmers))

    run_successfully, contigs = spades(spades_folder, threads, fastq_files, not_use_careful, max_memory,
                                       min_coverage_assembly, kmers, assembled_se_reads, spades_version=spades_version,
                                       estimated_coverage=estimated_coverage,
                                       spades_not_use_isolate=spades_not_use_isolate)

    if run_successfully:
        scaffolds = os.path.join(spades_folder, 'scaffolds.fasta')
        if keep_scaffolds:
            if os.path.isfile(scaffolds):
                shutil.copyfile(scaffolds, os.path.join(outdir, str('SPAdes_original_assembly.scaffolds.fasta')))
            else:
                print('The scaffolds file was not found!')

        if os.path.isfile(contigs):
            shutil.copyfile(contigs, os.path.join(outdir, str('SPAdes_original_assembly.contigs.fasta')))

            contigs_link = os.path.join(outdir, str(sample_name + '.contigs.fasta'))
            os.symlink(contigs, contigs_link)

            contigs = contigs_link

            min_contigs_length = define_minContigsLength(maximum_reads_length, min_contigs_length)

            sequence_dict = get_SPAdes_sequence_information(contigs)

            warnings, sequence_dict, filtered_sequences_sufix, spades_report_general = \
                decide_filter_parameters(sequence_dict, min_contigs_length, min_coverage_contigs,
                                         estimated_genome_size_mb, max_number_contigs)

            if filtered_sequences_sufix is not None:
                filtered_sequence_file = os.path.splitext(contigs)[0] + '.' + filtered_sequences_sufix + '.fasta'
                write_filtered_sequences_and_stats(sequence_dict, spades_report_general, contigs,
                                                   filtered_sequence_file, sample_name, False, save_excluded_contigs)
                contigs = filtered_sequence_file
            else:
                filtered_sequence_file = os.path.splitext(contigs)[0] + '.original.fasta'
                write_filtered_sequences_and_stats(sequence_dict, spades_report_general, contigs,
                                                   filtered_sequence_file, sample_name, True, False)
                contigs = filtered_sequence_file

            os.remove(contigs_link)
        else:
            run_successfully = False
            failing['sample'] = 'Assembly was not produced'
    else:
        failing['sample'] = 'Did not run'

    if not run_successfully:
        print(failing['sample'])
        pass_qc = False

    utils.removeDirectory(spades_folder)

    return run_successfully, pass_qc, failing, contigs, warnings
