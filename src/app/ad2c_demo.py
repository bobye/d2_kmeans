from ConfigParser import ConfigParser
import sys

config = ConfigParser(allow_no_value=True)
config.readfp(open(sys.argv[1]))



data_types_map={'EuclideanDist':'0', #default
                'Word2Vec':'7',
                'DenseHistogram':'5',
                'SparseHistogram':'12'}

options=''


options+=config.get('AD2c', 'executable')
options+=' -i ' + config.get('AD2c', 'input')

if config.get('AD2c', 'data') != '':
    options+=' --types ' + ','.join([ data_types_map[typename] for typename in config.get('AD2c', 'data').split(',')])
        
if config.get('AD2c', 'meta') != '':
    if (config.get('AD2c', 'phases') is None) or (config.get('AD2c', 'phases') is '1'):
        options+=' --metafile ' + config.get('AD2c', 'meta')
    else:
        print "Warning: have more than one phases, metafile ignored!"

if config.get('AD2c', 'phases') !='':
    options+=' -p ' + config.get('AD2c', 'phases')

if config.get('AD2c', 'dims') !='':
    options+=' -d ' + config.get('AD2c', 'dims')

if config.get('AD2c', 'supps') !='':
    options+=' -s ' + config.get('AD2c', 'supps')

if config.get('AD2c', 'clusters') != '':
    options+=' -k ' + config.get('AD2c', 'clusters')

if config.get('AD2c', 'output') != '':
    options+=' -o ' + config.get('AD2c', 'output')
        
if config.get('AD2c', 'cpus') == '1':
    comm=options
else:
    comm=options + ' -n ' + config.get('AD2c', 'size')
    comm+=' --prepare_batches ' + config.get('AD2c', 'cpus')
    comm+=' # split data cpu-wise\n'
    comm+='mpirun -n ' + config.get('AD2c', 'cpus') + ' ' + options
    comm+=' -n ' + str((int(config.get('AD2c', 'size'))-1) / int(config.get('AD2c', 'cpus')) + 1)
    comm+=' # run mpi job\n'

    
if config.get('AD2c', 'is_test_data') == 'Yes':
    if (config.get('AD2c', 'centroid') != ''):
        comm+=' --eval ' + config.get('AD2c', 'centroid')
    else:
        print "Warning: centroid file missing!"

print comm

