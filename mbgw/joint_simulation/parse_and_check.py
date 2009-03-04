# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from csv import reader
from numpy import ndarray, array, isscalar
import sys

def line_to_desc(line):
    while True:
        try:
            line.remove('')
        except:
            break
    return int(line[0]), line[1], line[2]

def translate(cls):
    if cls=='NULL':
        return type(None)
    elif cls=='scalar':
        return 'scalar'
    elif cls=='vector':
        return 'vector'
    elif cls=='matrix':
        return ndarray
    

def parse_tree(tree_reader):
    len, cls, key = line_to_desc(tree_reader.next())

    # Leaf case
    if cls != 'list':
        return translate(cls), key
    first_obj, first_key = parse_tree(tree_reader)

    # List case
    if first_key == 'NONE':
        out = [first_obj]
        for i in xrange(len-1):
            this_obj, this_key = parse_tree(tree_reader)
            out.append(this_obj)
        return out, key
    
    # Dict case    
    else:
        out = {first_key: first_obj}
        for i in xrange(len-1):
            this_obj, this_key = parse_tree(tree_reader)
            out[this_key] = this_obj
        return out, key

class MatchingError(TypeError):
    pass
            
def compare_tree(obj_tree, class_tree):
    
    # List case
    if isinstance(obj_tree, list):
        if class_tree == 'vector':
            if isinstance(obj_tree, list) or isinstance(obj_tree, ndarray):
                return obj_tree
            else:
                raise MatchingError, ': Classes do not match: should be %s but is %s'%(class_tree, type(obj_tree))
        if not isinstance(class_tree, list):
            raise MatchingError, ': Classes do not match: should be %s but is %s'%(class_tree, type(obj_tree))
        new_tree = []
        for i in xrange(len(obj_tree)):
            try:
                new_tree.append(compare_tree(obj_tree[i], class_tree[i]))
            except MatchingError:
                cls, inst, traceback = sys.exc_info()
                msg = '[%i]'%i + inst.message
                raise MatchingError, msg, traceback
        return new_tree
        
    # Dict case
    elif isinstance(obj_tree, dict):
        if not isinstance(class_tree, dict):
            raise MatchingError, ': Classes do not match: should be %s but is %s'%(class_tree, type(obj_tree))
        new_tree = {}
        for key in obj_tree.iterkeys():
            try:
                new_tree[key] = compare_tree(obj_tree[key], class_tree[key])
            except MatchingError:
                cls, inst, traceback = sys.exc_info()
                msg = "['%s']"%key + inst.message
                raise MatchingError, msg, traceback                            
        return new_tree
    
    # Scalar case
    if class_tree == 'scalar':
        if  isscalar(obj_tree):
            return obj_tree
    
    # Leaf case
    elif type(obj_tree) != class_tree:
        if isinstance(class_tree, list) and len(class_tree) == 1:
            if type(obj_tree) == class_tree[0]:
                return [obj_tree]
        raise MatchingError, ': Classes do not match: should be %s but is %s'%(class_tree, type(obj_tree))
    else:
        return obj_tree

if __name__ == '__main__':    
    class_tree, krap = parse_tree(reader(file('listSummary_preLoopObj_month1.txt'), delimiter=' '))
    obj_tree = {
    'KrigMatObj': {'PostMeanInterim.FULLlist': array(3),
                    'PostMeanInterim.LISTS': [None,
                                              array(3),
                                              array(3),
                                              array(3)],
                    'PostVar.FULLlist': array(3),
                    'PostVar.LISTS': [None,
                                      array(3),
                                      array(3),
                                      array(3)],
                    'cPtoP': array(3)},
     'OutMATlist': [array(3), 1],
     'footprintObj': {'footprintLUT': array(3),
                      'footprintXYT': array(3),
                      'modifiedFootprintTable': array(3),
                      'tseq': 1,
                      'xseq': 2.,
                      'yseq': 1}}

    class_tree['KrigMatObj']['PostVar.FULLlist'] = [ndarray]
    new_obj_tree = compare_tree(obj_tree,class_tree)