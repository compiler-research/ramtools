"""Usage: tools_perf.py generate [-n NUMBER] [--out OUTFILE] GENOMETABLE
          tools_perf.py convert bam [-iI] [--out OUTFOLDER] SAMFILE [BAMFILE]
          tools_perf.py convert ram [-iINrT] [-a ALG] [--out OUTFOLDER] SAMFILE ROOTFILE
          tools_perf.py index bam [-iI] [--out OUTFOLDER] FILE...
          tools_perf.py index ram [-iINS] [--out OUTFOLDER] FILE...
          tools_perf.py view bam FILE VIEWS [-iIP] [--out OUTFOLDER] [RANGE]
          tools_perf.py view ram FILE VIEWS [-ciINPsT] [--out OUTFOLDER] [RANGE] [--macro MACRO]
          tools_perf.py parsetreestats TTREEPERFSTATS...

Preprocessing and postprocessing for evaluating the performance of ramtools functions

Arguments:
  GENOMETABLE       CSV file with name of genome and ranges for each RNAME
  VIEWS             CSV file with genome, rname and region
  RANGE             Range in comma/dash separated value for the experiments to execute
  LOGFILE           Output of calling the tools_perf run on a set of files
  FILE              ROOTfile for the provided genome
  TTREEPERFSTATS    File with TTreePerfStats

Options:
  -h --help
  -a, --alg ALG         Compression algorithm of choice
  -c --cache            Enable TTreeCache
  -i --interactive      Interactive mode, prints commands before running them
  -I --io               Profile IO operations
  -m, --macro MACRO     Custom ramview macro to crossvalidate
  -n, --num NUMBER      Amount of records to generate
  -N --no-compile       Avoid compiling code
  -o, --out OUTFILE     File to save/append values, defaults to stdin
  -P, --parallel        Run scripts in parallel
  -r --no-index         Convert as raw without index
  -s --stats            Print TTreeStats to file
  -S --separate         Put index in separate file
  -T --no-split         Reduce Splitlevel for banches
"""
import os
import random
import sys
from docopt import docopt
import pandas as pd
import subprocess
from datetime import datetime

processes = []


def rangestr2list(s):
    return sum(((list(range(*[int(j) + k for k, j in enumerate(i.split('-'))]))
                if '-' in i else [int(i)]) for i in s.split(',')), [])


def clear_buffer_cache():
    code = subprocess.call(['sudo', 'sh', '-c', "echo 3 > /proc/sys/vm/drop_caches"])
    if code == 0:
        return
    else:
        print(code)
        sys.exit(1)


def sudo_reauth():
    subprocess.call('while true; do sudo -n true; sleep 60; kill -0 "$$" || exit; done 2>/dev/null &', shell=True)


def lauch_and_save_output(cmd, outfile, operation=None, interactive=False):
    if operation is None:
        operation = cmd

    print_timestamp()
    print("Executing {0}".format(operation))

    if interactive:
        manual_check(cmd)

    with open(outfile, 'w') as f:
        processes.append(subprocess.Popen(cmd, stdout=f))


def manual_check(cmd):
    print(" ".join(cmd))
    print("\nIs this the command you want yo issue [y/N]")
    while(True):
        choice = input().lower()

        if choice == 'y':
            break
        elif choice == 'n':
            sys.exit(1)
        else:
            print("Invalid Choice [y/N]")


def print_timestamp():
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print("[{0}] ".format(timestamp), end="")


def wait_for_all(processes):
    exit_codes = [p.wait() for p in processes]
    return exit_codes


def wrap_root_cmd(cmd):
    return ["root.exe", "-q", "-l", "-b"] + cmd


def wrap_time_cmd(cmd, time_logfile):
    return ["/usr/bin/time", "-v", "--output={0}".format(time_logfile)] + cmd


def wrap_io_cmd(cmd, io_logfile):
    return ['strace', '-o', io_logfile, '-TC'] + cmd


if __name__ == '__main__':
    arguments = docopt(__doc__)
    compilation_flag = '+' if not arguments['--no-compile'] else ''

    if arguments['generate']:

        N = int(arguments['--num']) if arguments['--num'] else 10

        outfile = arguments['--out'] if arguments['--out'] else None

        offset = 0
        if outfile is not None:
            if not os.path.isfile(outfile):
                with open(outfile, 'w') as f:
                    print(",rname,start,end", file=f)
            else:
                df = pd.read_csv(outfile)
                offset = df.index.max() + 1

        outfile = open(outfile, 'a') if outfile is not None else sys.stdout

        table = pd.read_csv(arguments['GENOMETABLE'])
        table = table[~table['RNAME'].str.startswith('GL')]
        table = table[~table['RNAME'].str.startswith('*')]

        for i in range(N):
            row = table.ix[random.choice(table.index)]
            rname = row['RNAME']
            a, b = random.randint(row['beginPOS'], row['endPOS']), random.randint(row['beginPOS'], row['endPOS'])
            a, b = min(a, b), max(a, b)
            print("{0},{1},{2},{3}".format(i+offset, rname, a, b), file=outfile)

    else:
        outfolder = arguments['--out'] if arguments['--out'] else '.'
        os.makedirs(outfolder, exist_ok=True)

        if arguments['convert']:
            if arguments['bam']:

                samfile = arguments['SAMFILE']
                bamfile = arguments['BAMFILE'] if arguments['BAMFILE'] else arguments['SAMFILE'].replace('.sam', '.bam')

                sam_basename = os.path.basename(samfile).split('.sam')[0]
                logfile = "samtobam_{0}".format(sam_basename)
                logfile = os.path.join(outfolder, logfile)

                convert_cmd = ['samtools', 'view', '-bS', '-o', bamfile, samfile]
                operation = "samtobam from {0} to {1}".format(samfile, bamfile)

            elif arguments['ram']:
                samtoram_macro = arguments['--macro'] if arguments['--macro'] else "samtoram.C"

                samfile = arguments['SAMFILE']
                rootfile = arguments['ROOTFILE']

                split = "false" if arguments['--no-split'] else "true"
                index = "false" if arguments['--no-index'] else "true"
                compression = arguments['--alg'] if arguments['--alg'] else "ROOT::kLZMA"

                split_name = 'split' if split == 'true' else 'nosplit'
                index_name = 'index' if index == 'true' else 'noindex'
                compression_name = compression.split('ROOT::k')[1].lower()
                sam_basename = os.path.basename(samfile).split('.sam')[0]

                logfile = "samtoram_{0}_{1}_{2}_{3}".format(sam_basename, compression_name, split_name, index_name)
                logfile = os.path.join(outfolder, logfile)

                convert_cmd = ['samtoram.C{5}("{0}", "{1}", {2}, {3}, {4})'.format(samfile, rootfile, index,
                                                                                   split, compression, compilation_flag)]
                convert_cmd = wrap_root_cmd(convert_cmd)
                operation = "samtoram from {0} to {1}".format(samfile, rootfile)

            convert_cmd = wrap_time_cmd(convert_cmd, logfile + '.perf')
            if arguments['--io']:
                convert_cmd = wrap_io_cmd(convert_cmd, logfile + '.io')

            lauch_and_save_output(convert_cmd, logfile + ".log", operation, interactive=arguments['--interactive'])
            wait_for_all(processes)

            if arguments['bam']:
                with open(logfile + ".log", 'w') as f:
                    print(subprocess.check_output(['stat', bamfile]).decode()[:-1], file=f)
                    print(subprocess.check_output(['stat', samfile]).decode(), file=f)

        if arguments['index']:
            sudo_reauth()

            for file in arguments['FILE']:

                clear_buffer_cache()

                file_basename = os.path.basename(file).split('.sam')[0]
                logfile = "index__{0}".format(file_basename)
                logfile = os.path.join(outfolder, logfile)

                if arguments['bam']:

                    logfile = logfile.replace('index__', 'bamindex__')

                    index_cmd = ['samtools', 'index', file, file + '.bai']
                    operation = "samtools index on {0}".format(file)

                elif arguments['ram']:

                    logfile = logfile.replace('index__', 'ramindex__')

                    index_cmd = ['ramindex.C{1}("{0}")'.format(file, compilation_flag)]

                    if arguments['--separate']:
                        index_cmd[-1] = index_cmd[-1][:-1] + ', false)'

                    index_cmd = wrap_root_cmd(index_cmd)

                    operation = "ramtools index {0}".format(file)

                index_cmd = wrap_time_cmd(index_cmd, logfile + '.perf')
                if arguments['--io']:
                    index_cmd = wrap_io_cmd(index_cmd, logfile + '.io')

                lauch_and_save_output(index_cmd, logfile + ".log", operation, interactive=arguments['--interactive'])

                wait_for_all(processes)

                if arguments['bam']:
                    with open(logfile + ".log", 'w') as f:
                        print(subprocess.check_output(['stat', file + '.bai']).decode(), file=f)

        elif arguments['view']:
            df = pd.read_csv(arguments['VIEWS'])

            if arguments['RANGE']:
                df = df.ix[rangestr2list(arguments['RANGE'])]

            sudo_reauth()
            clear_buffer_cache()

            for index, row in df.iterrows():

                region = "{0}:{1}-{2}".format(row['rname'], row['start'], row['end'])

                file = arguments['FILE'][0]
                logfile = "OP__{0}__{1}".format(os.path.basename(file), region)
                logfile = os.path.join(outfolder, logfile)

                if arguments['bam']:
                    logfile = logfile.replace("OP__", "bamview__")
                    cmd = ["samtools", "view", file, region]
                    operation = "samtools view on {0} {1}".format(file, region)

                elif arguments['ram']:
                    logfile = logfile.replace("OP__", "ramview__")

                    # Options
                    ramview_macro = arguments['--macro'] if arguments['--macro'] else "ramview.C"
                    cache = "false" if arguments['--cache'] else "true"

                    cmd = ['{2}{4}("{0}", "{1}", {3})'.format(file, region, ramview_macro, cache, compilation_flag)]
                    cmd = wrap_root_cmd(cmd)

                    if arguments['--stats']:
                        ttreeperf_file = logfile + '.root'
                        cmd[-1] = cmd[-1][:-1] + ', true, "{0}")'.format(ttreeperf_file)
                    operation = "ramtools view on {0} {1}".format(file, region)

                if not arguments['--io']:
                    cmd = wrap_time_cmd(cmd, logfile + '.perf')

                if arguments['--io']:
                    cmd = wrap_io_cmd(cmd, logfile + '.io')

                lauch_and_save_output(cmd, logfile + '.log', operation, interactive=arguments['--interactive'])

                if not arguments['--parallel']:
                    wait_for_all(processes)
                    clear_buffer_cache()

            wait_for_all(processes)

        elif arguments['parsetreestats']:
            for file in arguments['TTREEPERFSTATS']:
                if not file.endswith('.root'):
                    print("{0} does not have root extension, skiping...")
                    continue
                textfile = os.path.splitext(file)[0] + ".treeperf"
                imagefile = os.path.splitext(file)[0] + ".png"

                parsetreestats_cmd = ['parsetreestats.C+("{0}", true, "{1}")'.format(file, imagefile)]
                parsetreestats_cmd = wrap_root_cmd(parsetreestats_cmd)

                operation = "parsetreestats on {0}".format(file)
                lauch_and_save_output(parsetreestats_cmd, textfile, operation, interactive=arguments['--interactive'])

