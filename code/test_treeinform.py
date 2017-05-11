import os
import treeinform as ti
from ete3 import Tree
from collections import defaultdict

test3=Tree("((A@0:0.005, B@1:0.005):0.005,(A@2:0.005,B@3:0.007):0.005):1;")
# code will not work unless threshold is set above 0.33

def test_tree1():
    test1=Tree("(A@0:0.2, A@1:0.2):1;")

    ti.annotate_node(test1, 0.05)
    assert(test1.under==False)
    candidates=ti.prune(test1)
    assert(bool(candidates)==False)

    ti.annotate_node(test1, 0.41)
    assert(test1.under==True)
    candidates=ti.prune(test1)
    assert(candidates[('A', 0.4)]==set(['1', '0']))
    ti.output(candidates,'test.out1')
    assert(os.system('cat test.out1')=='0.4 1 0' or '0.4 0 1')

def test_tree2():
    test2=Tree("((B@0:0.2, A@1:0.2):0.1,A@2:0.4):1;")

    ti.annotate_node(test2, 0.45)
    assert(test2.under==False)
    candidates=ti.prune(test2)
    assert(candidates[('A',0.4)]==set(['1']))
    assert(candidates[('B',0.4)]==set(['0']))
    ti.output(candidates, 'test.out2')
    assert(os.system('cat test.out2')==0)

    ti.annotate_node(test2, 0.95)
    assert(test2.under==True)
    candidates=ti.prune(test2)
    assert(candidates[('A',0.9)]==set(['1','2']))
    ti.output(candidates, 'test.out2')
    assert(os.system('cat test.out2')=='0.9 1 2' or '0.9 2 1')

def main_test():
    test_tree="[&U] (((((A@0:0.01,A@1:0.01):0.01,B@3:0.01):0.01,A@2:0.01):0.01,B@4:0.01):0.8,(B@5:0.6,C@8:0.4)):1;"
    with open('test_tree', 'w') as f: f.write(test_tree)
    os.system('python treeinform.py -t test_tree -b 0.1 -o test.out')
    assert(os.system('cat test.out | wc -l')==2)
