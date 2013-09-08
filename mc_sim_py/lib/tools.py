'''

@author: Jakub Krajniak <jkrajniak@gmail.com>
'''

import os
import cPickle, gzip
import logging
import bisect

import numpy as np

def interaction_table_to_dict(file_name):
    """
    Return amino acids interaction table
    """
    return_dict = {}
    
    lines = [ l.strip("\n").strip().split(" ") for l in open(file_name, "r").readlines() ]
    aa_list = lines[0]
    for line in lines[1:]:
        
        return_dict[line[0]] = dict(zip(aa_list, map(float, line[1:len(line)-1])))
    
    return return_dict

def build_monomers(file_name):
    """
    Build monomers, based on interaction table
    """
    class_template = """
class %s(Monomer):
    name = "%s"
    interaction_table = %s
    \n\n
"""
    
    interaction_dict = interaction_table_to_dict(file_name)

    if os.path.exists(file_name):
        answer = raw_input('File extists, overwrite [Y/N] ? ')
        if answer == 'N':
            return False

    class_file = open('Monomer.py', 'w+')
    for key in interaction_dict.keys():
        cls_name = "Type%s" % key
        table = interaction_dict[key]

        string = class_template % (cls_name, key, str(table))
        class_file.writelines(string)

    class_file.close()

def create_dir(file_name):
    try:
        os.makedirs(os.path.dirname(file_name))
        return True
    except OSError, e:
        dir_path = os.path.dirname(file_name)
        if dir_path and not os.path.exists(os.path.dirname(file_name)):
            raise e
        elif not os.path.exists(file_name):
            raise e

def load_file(file_name, default=None):
    try:
        v = cPickle.load(gzip.open(os.path.realpath(file_name), "rb"))
        print "Open:",os.path.realpath(file_name), "Type:", type(v)
        return v
    except Exception, e:
        logging.debug(e)
        return default

def save_file(file_name, obj):
    return cPickle.dump(obj=obj, file=gzip.open(file_name, "wb", compresslevel=3), protocol=2)
    #return cPickle.dump(obj = obj, file = open(file_name, "wb"), protocol = 2)

def clean_dir(dir_name):
    
    for the_file in os.listdir(dir_name):
        file_path = os.path.join(dir_name, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception, e:
            print e

def array_diff(array_a, array_b):
    """
    Return value difference for two arrays
    """
    return map(lambda x,y: x - y, array_a, array_b)

def dot(vector_a, vector_b = None):
    """
    Dot products of two vectors
    """
    if vector_b == None:
        vector_b = (0,0,0)
    
    distance = array_diff(vector_a, vector_b)    
    return sum(map(lambda x: x**2, distance))

def index(list_object, item):
    """
    Locate the leftmost value exactly equal to x
    """
    i = bisect.bisect_left(list_object, item)
    if i != len(list_object) and list_object[i] == item:
        return i
    raise ValueError

def generate_one_xyz_file(path):
    files = os.listdir(path)
    jmol_files = filter(lambda x: "_" in x, files)
    jmol_files.sort(key = lambda x: int(x.split("_")[2].split(".")[0]))
    
def chunk(l, n):
    return (l[i:i+n] for i in xrange(0, len(l), n))

def get_column(data, column):
    return [ float(x.split(';')[column]) for x in data ]

def f_avg(column):
    return np.average(column)

def f_avg2(column):
    return np.average(map(lambda x: x**2, column))

def f_cv(e, e2, temp):
    ret = []
    for idx in range(len(temp)):
        val = (1/(temp[idx]**2)) * (e2[idx] - (e[idx]**2))
        ret.append(val)

    return ret

def f_cv2(data_e, N):
    ret = []
    for temp, e in data_e.iteritems():
        avg_e = np.average(data_e[temp])
        avg_e2 = avg_e*avg_e
        tmp = (np.average([ np.math.pow((data_e[temp][x]), 2) for x in range(len(data_e[temp])) ]) - avg_e2) / (N*(temp**2))
        ret.append( (temp, tmp))
    print ret
    ret.sort(key = lambda x: x[0] )
    print ret
    r = [ x[1] for x in ret ]
    return r

def avg_dict(input):
    """
    Return for given temp, list of steps (if more than one, do average)
    """
    ret = {}
    for temp in input.keys():
        step_list = input[temp]
        step_keys = step_list.keys()
        step_keys.sort()
        for step in step_keys:
            t = input[temp][step]
            if temp in ret:
                ret[temp].append(np.average(t))
            else:
                ret[temp] = [np.average(t)]
    return ret

def save_dict_to_csv(input_dict, csv_name):
    """
    Save internal dicts to csv files
    """
    f = open(csv_name, "w+")
    lines = []
    steps = input_dict.keys()
    steps.sort()
    for step in steps:
        lines.append("%d;%s\n" % (step, ";".join(map(str, input_dict[step]))))

    f.writelines(lines)


def get_list(csv, delimiter=';', data_type=float):
    l = map(data_type, csv.split(delimiter))
    return l
