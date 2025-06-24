import os
import matplotlib.pyplot as plt
import numpy as np
import re
from glob import glob
import pandas as pd

columns = ['file', 'genome', 'method', 'index', 'cache', 'alg', 'split', 'memlim', 'region', 'usertime', 'systemtime', 'cpu_usage', 'memory', 'filesize']


def load_perf_file(file, method=""):
    if method == "":
        method = os.path.basename(file).split('_')[0]

    with open(file, 'r') as f:
        s = f.read()

    genome = os.path.basename(file).split('__')[1]
    region = os.path.basename(file).split('__')[2].split('.perf')[0]

    lines = s.split('\n')
    usertime = float(lines[1].split(': ')[1])
    systemtime = float(lines[2].split(': ')[1])
    cpu_usage = float(lines[3].split(': ')[1].split('%')[0])
    memory = int(lines[9].split(': ')[1].split('%')[0])

    filesize = os.stat(file.replace('.perf', '.log')).st_size

    if '_ram' in file:
        method = 'ramtools'
        if '_lzma' in file:
            alg = 'lzma'
        elif '_zlib' in file:
            alg = 'zlib'
        else:
            alg = ''

        split = '_nosplit' not in file
        index = '_index' in file
        cache = '_cache' in file
    elif '_sam' in file:
        method = 'samtools'
        alg = split = index = cache = None
    else:
        raise ValueError(method)

    memlim = 'memlim' in file

    basegenome = os.path.splitext(os.path.basename(genome))[0]
    if '_' in basegenome:
        basegenome = basegenome.split('_')[0]

    return [genome, basegenome, method, index, cache, alg, split, memlim, region, usertime, systemtime, cpu_usage, memory, filesize]


def load_perf_folder(folder):
    perfs = [load_perf_file(f) for f in glob('{0}/*.perf'.format(folder))]

    df = pd.DataFrame(data=perfs, columns=columns)
    df['totaltime'] = df['usertime'] + df['systemtime']

    df['Speed_MBps'] = (df['filesize']/1024**2) / df['totaltime']

    df['mem_MB'] = df['memory']/(1024)
    df['size_MB'] = (df['filesize']/(1024**2)).round(2)

    df = df.sort_values(['genome', 'method', 'index', 'cache', 'alg', 'split', 'memlim', 'region'])

    return df


def load_perf_superfolders(folders):
    df = pd.DataFrame(columns=columns)
    for folder in folders:
        df = pd.concat([df]+[load_perf_folder(subfolder) for subfolder in glob('{0}/*'.format(folder))], ignore_index=True)

    df = df.sort_values(['genome', 'method', 'index', 'cache', 'alg', 'split', 'memlim'])
    return df


def get_metric(df, column, regions):
    return np.array([df[df['region'] == r][column].values for r in regions])


def compare_metrics(df, methods, regions, column, save=False, relative=None, log=False):
    plt.figure(figsize=(25, 6))

    dfs = [df[df['method'] == m] for m in methods]
    metrics = [get_metric(d, column, regions) for d in dfs]

    x = np.arange(len(regions))
    N = len(methods) + 1

    for i, (m, method) in enumerate(zip(metrics, methods)):
        if relative is None:
            plt.bar(x+1/N*i, m, width=1/N, label=method)
        else:
            plt.bar(x+1/N*i, m/metrics[relative], width=1/N, label=method)
        # print(method, m)

    plt.xticks(x, regions, rotation=15)
    plt.title(column, fontsize=25)
    plt.legend(fontsize=16)

    if log:
        plt.yscale('log')

    if save:
        plt.savefig("images/{0}.png".format(column), format='png')


def load_samtoram_perf(file, method=''):

    with open(file, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    if "Command exited" in lines[0]:
        lines = lines[1:]

    if method == '':
        method = os.path.basename(lines[0].split(', ')[1]).strip(' "')

    usertime = float(lines[1].split(': ')[1])
    systemtime = float(lines[2].split(': ')[1])
    cpu_usage = float(lines[3].split(': ')[1].split('%')[0])
    memory = int(lines[9].split(': ')[1].split('%')[0])

    logfile = file.replace('.perf', '.log')

    with open(logfile, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    filesize = int(lines[4].split(' = ')[-1].strip(' *'))
    compression = float(lines[5].split(' = ')[-1].strip(' *'))

    file = method.split('_')[0]
    method = "ramtools_" + method[len(file):].split('.root')[0]

    return [file, method, usertime, systemtime, cpu_usage, memory, filesize, compression]


def load_samtobam_perf(file, method=''):

    with open(file, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    if "Command exited" in lines[0]:
        lines = lines[1:]

    if method == '':
        method = lines[0].split(' ')[-1].strip(' "').replace('.sam', '.bam')

    usertime = float(lines[1].split(': ')[1])
    systemtime = float(lines[2].split(': ')[1])
    cpu_usage = float(lines[3].split(': ')[1].split('%')[0])
    memory = int(lines[9].split(': ')[1].split('%')[0])

    logfile = file.replace('.perf', '.log')

    with open(logfile, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    filesize = int(lines[1].split('\t')[0].split(': ')[1])
    filesize_org = int(lines[10].split('\t')[0].split(': ')[1])
    compression = filesize_org / filesize

    file = os.path.basename(method).split('.bam')[0]

    return [file, 'samtools', usertime, systemtime, cpu_usage, memory, filesize, compression]


def load_bamindex_perf(file, method=''):

    with open(file, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    if "Command exited" in lines[0]:
        lines = lines[1:]

    if method == '':
        method = os.path.basename(lines[0].split(' ')[-2])

    usertime = float(lines[1].split(': ')[1])
    systemtime = float(lines[2].split(': ')[1])
    cpu_usage = float(lines[3].split(': ')[1].split('%')[0])
    memory = int(lines[9].split(': ')[1].split('%')[0])

    logfile = file.replace('.perf', '.log')

    with open(logfile, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    filesize = int(lines[1].split('\t')[0].split(': ')[1])

    file = method.split('.bam')[0]

    return [file, 'samtools', usertime, systemtime, cpu_usage, memory, filesize]


def load_ramindex_perf(file, method=''):

    with open(file, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    if "Command exited" in lines[0]:
        lines = lines[1:]

    if method == '':
        method = os.path.basename(lines[0].split('"')[2])

    usertime = float(lines[1].split(': ')[1])
    systemtime = float(lines[2].split(': ')[1])
    cpu_usage = float(lines[3].split(': ')[1].split('%')[0])
    memory = int(lines[9].split(': ')[1].split('%')[0])

    logfile = file.replace('.perf', '.log')

    with open(logfile, 'r') as f:
        s = f.read()
    lines = s.split('\n')

    for line in lines:
        if line.startswith('Size'):
            filesize = int(line.split('\t')[0].split(': ')[1])

    file = method.split('_')[0]
    method = "ramtools_" + method[len(file):].split('.root')[0]

    return [file, method, usertime, systemtime, cpu_usage, memory, filesize]


def load_samtoram_folder(folder):
    perfs = []
    for file in glob('{0}/samtoram*.perf'.format(folder)):
        print(file)
        perfs += [load_samtoram_perf(file)]

    for file in glob('{0}/samtobam*.perf'.format(folder)):
        print(file)
        perfs += [load_samtobam_perf(file)]

    columns = ['file', 'method', 'usertime', 'systemtime', 'cpu_usage', 'memory', 'filesize', 'compression']
    df = pd.DataFrame(data=perfs, columns=columns)
    df['size_GB'] = df['filesize']/(1024**3)
    return df


def load_index_folder(folder):
    perfs = []
    for file in glob('{0}/ramindex*.perf'.format(folder)):
        print(file)
        perfs += [load_ramindex_perf(file)]

    for file in glob('{0}/bamindex*.perf'.format(folder)):
        print(file)
        perfs += [load_bamindex_perf(file)]

    columns = ['file', 'method', 'usertime', 'systemtime', 'cpu_usage', 'memory', 'filesize']
    df = pd.DataFrame(data=perfs, columns=columns)
    df['size_GB'] = df['filesize']/(1024**3)
    return df
