#!/usr/bin/python

import os, sys, re, time
import subprocess
import getopt

def simpe_cmd_output(cmdargs):
    p = subprocess.Popen(cmdargs, stdout=subprocess.PIPE,)
    ret = []
    while True:
        l = p.stdout.readline()
        if not l:
            break
        ret.append(l)
        print l,

    while p.poll() == None:
        time.sleep(0.1)

    return ret

def simple_runcmd(f,args):
    f.write('$ %s\n' % (' '.join(args)))
    f.write(simpe_cmd_output(args)[0])
    f.write('\n')

def checktimestamp_sitespecific(f):
    o = simpe_cmd_output(['ssh', 'knight', 'date', '+"%s %N"'])
    ts = time.time()

    tmp = o[0].split()
    ts2 = float(tmp[0]) + float(tmp[1])*1e-9

    f.write('[checking time]\n')
    simple_runcmd(f, ['date'])
    f.write( 'ts on this machine: %f\n' % ts )
    f.write( 'ts on knight: %f\n' % ts2 )
    f.write( 'ts diff: %f\n' % abs(ts2-ts))
    f.write('\n')

    if abs(ts2-ts) > 0.1:
        print 'Please re-synch time!'
        sys.exit(1)


def output_env(f, platform='Altera'):
    f.write('platform: %s\n' % platform)

    simple_runcmd(f, ['hostname', '-f'])

    if platform == 'Altra':
        simple_runcmd(f, ['aocl', 'version'])

    if platform == 'Intel':
        simple_runcmd(f, ['lscpu'])
        simple_runcmd(f, ['numactl', '-H'])
        simple_runcmd(f, ['cat', '/opt/intel/opencl-1.2-6.4.0.25/doc/version.txt'])

    simple_runcmd(f, ['git', 'rev-parse', 'HEAD'])

    simple_runcmd(f, ['cat', '/proc/version'])

    simple_runcmd(f, ['cat', '/proc/loadavg'])



def bench(outputdir, benchtype, globalsize, platform):
    btable = {}
    btable[0] = ['XSBench', 'aocx']
    btable[1] = ['XSBench', 'aocx'] 
    # //

    label = 't%d-g%d' %  (benchtype, globalsize)
    outputfn = '%s/res-ocl-xsbench-%s.txt' % (outputdir, label)
#    etracelogfn = '%s/etrace-ocl-xsbench-%s.txt' % (outputdir, label)

    // switch benchtype: cp btable[benchtype[1]]  .

    args  = []
#    args += ['echo']
#    args += ['./etrace/etrace', '-i', '0.1', '-o', etracelogfn]
    if platform == 'Altera':
        args += ['taskset', '1']
    if platform == 'Intel':
        args += ['numactl', '-m', '0',  '--cpunodebind=0']
    if platform == 'Intel1':
        args += ['numactl', '-m', '1',  '--cpunodebind=1']
    args += [ btable[benchtype][0] ]
    args += ['-t', '%d' % gs]

    print 'Command:', ' '.join(args)
    try:
        f = open(outputfn, 'w')
    except:
        print 'Failed to write to', outputfn
        sys.exit(1)


    p = subprocess.Popen(args, stdout=subprocess.PIPE,)

    while True:
        l = p.stdout.readline()
        if not l:
            break
        f.write(l)
        print l,

    while p.poll() == None:
        time.sleep(1)

    s = 'RC=%d' % p.poll()
    f.write( s + '\n')

    f.close()



def usage():
    print 'Usage:', sys.argv[0]
    print ''
    print '-g globalsize'
    print '-t benchtype'
    print '   0: base'
    print '   1: fma'
    print '   2: cc'
    print '   3: optbs'
    print '   4: vec'

if __name__ == '__main__':

    resdir = 'results'
    btype = 0 # BASE
    gsizes = [15000000] # no of look
    platform = 'Altera'


    try:
        opts, args = getopt.getopt(sys.argv[1:], 'ht:g:', [])
    except getopt.GetoptError, err:
        print err
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ('-h'):
            usage()
            sys.exit(0)
        elif o in ('-t'):
            btype = int(a)
        elif o in ('-g'):
            gsizes = map(lambda x: int(x), a.split(','))
        elif o in ('-d'):
            resdir = a
        elif o in ('-p'):
            platform = a

    label = 't%d-s%s' %  (btype, bufsizestr)

    envfn = '%s/env-%s.txt' % (resdir, label)
    try:
        f = open(envfn, 'w')
    except:
        print 'Failed to open:', envfn
        sys.exit(1)

    output_env(f, platform)

    if  platform == "Altera" :
        checktimestamp_sitespecific(f)

    f.close()

    #
    #
    ts1 = time.time()
    time.sleep(5)

    for gs in gsizes:
        bench(resdir, btype, gs, platform)

    time.sleep(5)
    ts2 = time.time()

    # site-specific
    # read WT310 data from db
    if  platform == "Altera" :
        wt310fn = '%s/wt310-duteros-ocl-iabench-%s.txt' % (resdir, label)

        dbcmd = '~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py'
        simpe_cmd_output(['ssh', 'knight', dbcmd, '%lf' % ts1, '%lf' % ts2, '/tmp/wt310.txt'])
        simpe_cmd_output(['scp', 'knight:/tmp/wt310.txt', wt310fn])

    sys.exit(0)
