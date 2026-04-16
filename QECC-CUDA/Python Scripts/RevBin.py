# def replace_multiple_strings(file_path, replacements):
#     """
#     Replaces multiple target strings with their corresponding replacement strings in a file.
#     :param file_path: Path to the .txt file
#     :param replacements: A dictionary where keys are target strings and values are the replacement strings
#     """
#     try:
#         # Step 1: Read the file into memory
#         with open(file_path, 'r') as file:
#             file_data = file.read()

#         # Step 2: Make all replacements in memory
#         for target_string, replacement_string in replacements.items():
#             file_data = file_data.replace(target_string, replacement_string)

#         # Step 3: Write the updated content back to the file
#         with open(file_path, 'w') as file:
#             file.write(file_data)

#         print(f"Replaced strings in {file_path}")

#     except FileNotFoundError:
#         print(f"Error: The file {file_path} was not found.")
#     except Exception as e:
#         print(f"An error occurred: {e}")

# replacements = {}

# for i in range(1, 64, 1):
#     bin_number = format(i, '06b')

#     rev_number = bin_number[: : -1]
    
#     replacements['(' + str(i) + ')'] = '(' + str(int(rev_number, 2)) + ')'
#     replacements['  '.join(bin_number)] = '  '.join(rev_number)

# # Gather all keys that are values of other key-value pairs
# keys_to_remove = {v for v in replacements.values() if v in replacements}

# # Remove those keys from the dictionary
# for key in keys_to_remove:
#     del replacements[key]

# print(replacements)

# # Usage example
# file_path = 'example.txt'

# replace_multiple_strings(file_path, replacements)

# %%
for i in range(1, 64, 1):
    print("Old: " + str(i))
    bin_number = format(i, '06b')

    rev_number = bin_number[: : -1]

    print("New: " + str(int(rev_number, 2)))

    print('  '.join(rev_number))


# %%
