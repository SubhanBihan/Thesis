# %%
# Imports and function defs. MUST be run beforehand.
import itertools
import cupy as cp
import csv
from functools import reduce

I = cp.array([[1, 0], [0, 1]])
X = cp.array([[0, 1], [1, 0]])
Z = cp.array([[1, 0], [0, -1]])
Y = cp.array([[0, -1j], [1j, 0]])

################## VARIABLES #################
# n of [[n, k, d]], and desired t
# S1 is leftmost stabiliser (to match ancilla number in paper)

# # [5, 1, 3]
# n = 5
# t = 1

# # Stabilisers
# S1 = reduce(cp.kron, [X, Z, Z, X, I])
# S2 = reduce(cp.kron, [I, X, Z, Z, X])
# S3 = reduce(cp.kron, [X, I, X, Z, Z])
# S4 = reduce(cp.kron, [Z, X, I, X, Z])
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

# Stabilisers
S1 = reduce(cp.kron, [I, I, I, X, X, X, X])
S2 = reduce(cp.kron, [I, X, X, I, I, X, X])
S3 = reduce(cp.kron, [X, I, X, I, X, I, X])
S4 = reduce(cp.kron, [I, I, I, Z, Z, Z, Z])
S5 = reduce(cp.kron, [I, Z, Z, I, I, Z, Z])
S6 = reduce(cp.kron, [Z, I, Z, I, Z, I, Z])
S_all = [S1, S2, S3, S4, S5, S6]
S1_str = "IIIXXXX"
S2_str = "IXXIIXX"
S3_str = "XIXIXIX"
S4_str = "IIIZZZZ"
S5_str = "IZZIIZZ"
S6_str = "ZIZIZIZ"
S_all_str = [S1_str, S2_str, S3_str, S4_str, S5_str, S6_str]

##############################################

mpauli = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
mp_err = {'X': X, 'Y': Y, 'Z': Z}

def get_permutations(s):
    return set(''.join(p) for p in itertools.permutations(s))

def t_errors_set(t, numqubits):
    t_errors = set()
    
    # permutations = set()
    permutations = set(''.join(p) for p in itertools.product(mp_err.keys(), repeat = t))
        
    for permutation in permutations:
        t_errors.update(get_permutations(f'{permutation}{'I' * (numqubits - t)}'))
    
    # for key in mp_err.keys():
    #     t_errors.update(get_permutations(f'{key}{'I' * (numqubits - t)}'))
        
    return t_errors

def errors_dict(err_str_set):
    t_errors = dict()
    
    for err_str in err_str_set:
        err_mat = 1
        for pauli_mat in err_str:
            err_mat = cp.kron(err_mat, mpauli[pauli_mat])
        
        t_errors[err_str] = err_mat
    
    return t_errors

def ancilla_measure(S, E):
    return '0' if cp.array_equal(cp.matmul(S, E), cp.matmul(E, S)) else '1'

anc = len(S_all)

# %%
####################################### Regular correction #########################################

Corr_list = []
indices = []


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
for i in range(1, 2**anc):
    if i not in indices:
        Corr_list.append([i, f'{'I' * n}'])

Sorted_Corr_list = sorted(Corr_list, key = lambda x: x[0])

Sorted_Corrs = [sublist[1] for sublist in Sorted_Corr_list]

# For .txt
with open('corrections.txt', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[word] for word in Sorted_Corrs])

print("TXT file has been generated.")

# For .csv
with open('corrections.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Sorted_Corr_list)
    # writer.writerows(Sorted_Corrs)

print("CSV file has been generated.")

# %%
####################################### Erasure correction #########################################

# def t_errors_set_D1(t_errs, numqubits):
#     t_errors_D1 = t_errors_set(t_errs, numqubits)
#     t_errors_D1 = {err for err in t_errors_D1 if all(err[j] == 'I' for j in range(2, numqubits))}

#     return t_errors_D1

def t_errors_set_D1(t, numqubits):
    t_errors_D1 = t_errors_set(t, numqubits)
    t_errors_D1 = {err for err in t_errors_D1 if err[0] != 'I'}

    return t_errors_D1

Corr_list = []
indices = []


for i in range(1, t + 1):
    i_errors_dict = errors_dict(t_errors_set_D1(i, n))
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
for i in range(1, 2**anc):
    if i not in indices:
        Corr_list.append([i, f'{'I' * n}'])

Sorted_Corr_list = sorted(Corr_list, key = lambda x: x[0])

Sorted_Corrs = [sublist[1] for sublist in Sorted_Corr_list]

# For .txt
with open('corrections_erasure.txt', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows([[word] for word in Sorted_Corrs])

print("TXT file has been generated.")

# For .csv
with open('corrections_erasure.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Sorted_Corr_list)
    # writer.writerows(Sorted_Corrs)

print("CSV file has been generated.")

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
    for i in range(1, 2**anc):
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
