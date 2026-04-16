# %%
base_str = "FastPauli::CArr<NUMQUBITS> cArr_"

for i in range(1, 7):
    for j in range(0, pow(2, i)):                # Will automatically exclude 2^i since range is endpoint exclusive. For same reason we use 7 even though 6 ancillas
        print(base_str + format(j, f'0{i}b') + ';')


# %%
# Python script to generate Stabiliser code for n (here n = 6) stages
NUMQUBITS = 'NUMQUBITS'
ket_list = ['ketH1P', 'ketH1N']  # ket variables for positive and negative states

def generate_binary_strings(n):
    """Generate all binary strings of length n."""
    return [bin(i)[2:].zfill(n) for i in range(2**n)]

for stage in range(1, 7):  # stages 1 to 6
    print(f"// =====({stage})=====")
    
    # First block (always use cArr_0 and cArr_1 in the first stage)
    if stage == 1:
        print(f"FastPauli::kron(cArr8.data(), cArrH1.data(), cArrCodewordtmp.data(), 1, {NUMQUBITS});")
        print(f"FiveQubitCode::S{stage}Sequential(cArr8.data(), err, errTeleport);")
        print(f"FiveQubitCode::applyAncillaKet(cArr_0.data(), 1, {NUMQUBITS}, {ket_list[0]}.data(), cArr8.data());")
        print(f"FiveQubitCode::applyAncillaKet(cArr_1.data(), 1, {NUMQUBITS}, {ket_list[1]}.data(), cArr8.data());")
    
    # Subsequent stages
    else:
        # Generate all combinations of binary strings of length (stage - 1)
        binary_strings = generate_binary_strings(stage - 1)
        for binary_str in binary_strings:
            print(f"FastPauli::kron(cArr8.data(), cArrH1.data(), cArr_{binary_str}.data(), 1, {NUMQUBITS});")
            print(f"FiveQubitCode::S{stage}Sequential(cArr8.data(), errTeleport);")
            print(f"FiveQubitCode::applyAncillaKet(cArr_{binary_str}0.data(), 1, {NUMQUBITS}, {ket_list[0]}.data(), cArr8.data());")
            print(f"FiveQubitCode::applyAncillaKet(cArr_{binary_str}1.data(), 1, {NUMQUBITS}, {ket_list[1]}.data(), cArr8.data());")
    
    print("\n")

# %%
for i in range(0, 64):
    print(f"cArrAggregate[{i}] = cArr_{i:06b};")
