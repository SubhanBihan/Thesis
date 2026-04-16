# %%
# Imports and function defs. MUST be run beforehand.
# import multiprocessing
# import cProfile
import itertools
import numpy as np
import csv
from functools import reduce
from collections import defaultdict
# from numba import jit
# import concurrent.futures
# import os
# import logging
# import traceback
# logging.basicConfig(level=logging.DEBUG)
# logger = logging.getLogger(__name__)


I = np.array([[1, 0], [0, 1]])
X = np.array([[0, 1], [1, 0]])
Z = np.array([[1, 0], [0, -1]])
Y = np.array([[0, -1j], [1j, 0]])

################## VARIABLES #################
# n of [[n, k, d]], and desired t
# S1 is leftmost stabiliser (to match ancilla number in paper)

# # [5, 1, 3]
# n = 5
# t = 1
# # FOR ERASURE CORRECTION
# max_errs = 2
# valid_errs = 1

# # Stabilisers
# S1 = reduce(np.kron, [X, Z, Z, X, I])
# S2 = reduce(np.kron, [I, X, Z, Z, X])
# S3 = reduce(np.kron, [X, I, X, Z, Z])
# S4 = reduce(np.kron, [Z, X, I, X, Z])
# S_all = [S1, S2, S3, S4]
# # S_all = [S4, S3, S2, S1]
# S1_str = "XZZXI"
# S2_str = "IXZZX"
# S3_str = "XIXZZ"
# S4_str = "ZXIXZ"
# S_all_str = [S1_str, S2_str, S3_str, S4_str]

# [7, 1, 3]
n = 7
t = 2
# FOR ERASURE CORRECTION
max_errs = 3
valid_errs = 1

# Stabilisers
S1 = reduce(np.kron, [I, I, I, X, X, X, X])
S2 = reduce(np.kron, [I, X, X, I, I, X, X])
S3 = reduce(np.kron, [X, I, X, I, X, I, X])
S4 = reduce(np.kron, [I, I, I, Z, Z, Z, Z])
S5 = reduce(np.kron, [I, Z, Z, I, I, Z, Z])
S6 = reduce(np.kron, [Z, I, Z, I, Z, I, Z])
S_all = [S1, S2, S3, S4, S5, S6]
S1_str = "IIIXXXX"
S2_str = "IXXIIXX"
S3_str = "XIXIXIX"
S4_str = "IIIZZZZ"
S5_str = "IZZIIZZ"
S6_str = "ZIZIZIZ"
S_all_str = [S1_str, S2_str, S3_str, S4_str, S5_str, S6_str]

# # [8, 3, 3]
# n = 8
# t = 2
# # FOR ERASURE CORRECTION
# max_errs = 3
# valid_errs = 1

# # Stabilisers
# S1 = reduce(np.kron, [X, I, X, I, Z, Y, Z, Y])
# S2 = reduce(np.kron, [Z, Y, I, X, Z, Y, I, X])
# S3 = reduce(np.kron, [I, X, I, X, Y, Z, Y, Z])
# S4 = reduce(np.kron, [I, X, Z, Y, I, X, Z, Y])
# S5 = reduce(np.kron, [I, Y, X, Z, X, Z, I, Y])
# S_all = [S1, S2, S3, S4, S5]
# S1_str = "XIXIZYZY"
# S2_str = "ZYIXZYIX"
# S3_str = "IXIXYZYZ"
# S4_str = "IXZYIXZY"
# S5_str = "IYXZXZIY"
# S_all_str = [S1_str, S2_str, S3_str, S4_str, S5_str]

##############################################

anc = len(S_all)
dim_anc = 2 ** anc

mpauli = {'X': X, 'Y': Y, 'Z': Z, 'I': I}
mp_err = {'X': X, 'Y': Y, 'Z': Z}

def get_permutations(s):
    return set(''.join(p) for p in itertools.permutations(s))

# Returns a set of error strings with t errors (in each)
def t_errors_set(t, numqubits):
    if t == 0:
        return {f'{'I' * numqubits}'}
    
    t_errors = set()
    
    permutations = set(''.join(p) for p in itertools.product(mp_err.keys(), repeat = t))
        
    for permutation in permutations:
        t_errors.update(get_permutations(f'{permutation}{'I' * (numqubits - t)}'))
    
    # for key in mp_err.keys():
    #     t_errors.update(get_permutations(f'{key}{'I' * (numqubits - t)}'))
        
    return t_errors

# Returns a dictionary corresponding to given error strings set, with key (str) : value (matrix)
def errors_dict(err_str_set):
    t_errors = dict()
    
    for err_str in err_str_set:
        err_mat = 1
        for pauli_mat in err_str:
            err_mat = np.kron(err_mat, mpauli[pauli_mat])
        
        t_errors[err_str] = err_mat
    
    return t_errors

# Gives ancilla qubit value provided a stabiliser and error matrix
def ancilla_measure(S, E):
    return '0' if np.array_equal(np.matmul(S, E), np.matmul(E, S)) else '1'


# %%
####################################### Regular correction #########################################

Corr_list = [None for _ in range(dim_anc)]


for i in range(1, t + 1):
    i_errors_dict = errors_dict(t_errors_set(i, n))
    for err_str, err_mat in i_errors_dict.items():
        list_mem = []
        ancilla_str = ""
        
        for S in S_all:
            ancilla_str += ancilla_measure(S, err_mat)
            
        index = int(ancilla_str, 2)
        if not Corr_list[index]:
            Corr_list[index] = err_str

# Initialise IIII... at remaining indices
for i in range(1, dim_anc):
    if not Corr_list[i]:
        Corr_list[i] = f'{'I' * n}'

Corr_list.pop(0)
Sorted_Corr_list = list(enumerate(Corr_list, start = 1))


# For .txt
with open('corrections.txt', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[word] for word in Corr_list])

print("TXT file has been generated.")

# For .csv
with open('corrections.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Sorted_Corr_list)
    # writer.writerows(Sorted_Corrs)

print("CSV file has been generated.")

# %%
####################################### Erasure correction #########################################

def Dk_err(t, numqubits, k, err):
    if err == 'I':
        t -= 1
    
    t_errs_set = t_errors_set(t, numqubits)
    
    return {k_err for k_err in t_errs_set if k_err[k] == err}


full_list = []

# def code_to_profile() -> None:
for k in range(n):
    # table = [[] for _ in range(dim_anc)]
    # Dk_normal_errs = []
    table = [None for _ in range(dim_anc)]
    counts = defaultdict(int)
    
    while any(elem is None for elem in table):
        for i in range(1, max_errs + 1):
            i_errors_dict = dict()
            for err in mpauli.keys():
                i_errors_dict.update(errors_dict(Dk_err(i, n, k, err)))
            
            # if i <= valid_errs:
            #     Dk_normal_errs.extend(list(i_errors_dict.keys()))
            
            for err_str, err_mat in i_errors_dict.items():
                ancilla_str = ""
                
                for S in S_all:
                    ancilla_str += ancilla_measure(S, err_mat)
                    
                index = int(ancilla_str, 2)
                
                if (not table[index]) and counts[err_str[k]] < dim_anc / 4:
                    table[index] = err_str
                    counts[err_str[k]] += 1
                
    
    table.pop(0)
    full_list.extend(list(enumerate(table, start = 1)))
    

Sorted_Corrs = [sublist[1] for sublist in full_list]

# For .txt
with open('corrections_erasure.txt', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[word] for word in Sorted_Corrs])

print("TXT file has been generated.")

# For .csv
with open('corrections_erasure.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(full_list)
    # writer.writerows(Sorted_Corrs)

print("CSV file has been generated.")

# cProfile.run("code_to_profile()", sort="tottime")

# %%
####################################### Flag correction #########################################

def second_perr_index(s: str) -> int:
    # Find the index of the 2nd character that is not 'I'
    count = 0
    for i, char in enumerate(s):
        if char != 'I':
            count += 1
        if count == 2:
            return i
    return -1  # Return -1 if there are fewer than 2 non-'I' characters

def second_last_perr_index(s: str) -> int:
    # Find the index of the 2nd-last character that is not 'I'
    count = 0
    for i in range(len(s) - 1, -1, -1):
        if s[i] != 'I':
            count += 1
        if count == 2:
            return i
    return -1  # Return -1 if there are fewer than 2 non-'I' characters

def perr_count(s: str, start_index: int) -> int:
    count = 0
    
    for i in range(start_index, len(s)):
        if s[i] != 'I':
            count += 1

    return count

def prop_errors(S_str: str):
    start = second_perr_index(S_str)
    stop = second_last_perr_index(S_str)
    prop_errors_set = set()
    total_perr = perr_count(S_str, 0)
    
    for i in range(start, stop + 1):
        base_str = S_str
        
        # Fixing minimal-weight base error
        if S_str[i] != 'I':
            if perr_count(S_str, i) > total_perr / 2:
                base_str = base_str[: i] + f'{'I' * (len(S_str) - i)}'
            else:
                base_str = f'{'I' * i}' + base_str[i :]
        
        prop_errors_set.update({base_str[: i] + char + base_str[i + 1:] for char in mpauli.keys()})

    return prop_errors_set


def flag_table(S_str: str):
    Corr_list = []
    indices = []
    
    prop_errors_dict = errors_dict(prop_errors(S_str))
    for err_str, err_mat in prop_errors_dict.items():
        list_mem = []
        ancilla_str = ""
        
        for S in S_all:
            ancilla_str += ancilla_measure(S, err_mat)
            
        index = int(ancilla_str, 2)
        if index in indices or index == 0:
            continue
        
        indices.append(index)
        list_mem.append(index)
        list_mem.append(err_str)
        
        Corr_list.append(list_mem)

    for i in range(1, t + 1):
        i_errors_dict = errors_dict(t_errors_set(i, n))
        for err_str, err_mat in i_errors_dict.items():
            list_mem = []
            ancilla_str = ""
            
            for S in S_all:
                ancilla_str += ancilla_measure(S, err_mat)
                
            index = int(ancilla_str, 2)
            if index in indices:
                continue
            
            indices.append(index)
            list_mem.append(index)
            list_mem.append(err_str)
            
            Corr_list.append(list_mem)
    
    # Initialise IIII... at remaining indices
    for i in range(1, dim_anc):
        if i not in indices:
            Corr_list.append([i, f'{'I' * n}'])
    
    Sorted_Corr_list = sorted(Corr_list, key = lambda x: x[0])
    
    return Sorted_Corr_list

Complete_corrs_list = []
           
for S in S_all_str:
    Complete_corrs_list.extend(flag_table(S))

Sorted_Corrs = [sublist[1] for sublist in Complete_corrs_list]

# For .txt
with open('corrections_flag.txt', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[word] for word in Sorted_Corrs])

print("TXT file has been generated.")

# For .csv
with open('corrections_flag.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Complete_corrs_list)
    # writer.writerows(Complete_corrs_list)

print("CSV file has been generated.")
# %%
