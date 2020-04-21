import os
# # for dirpath, dirname, filename in os.walk()
# for path in snakemake.input:
#     for dirpath, dirnames, filenames in os.walk(path):
#         if os.path.basename(dirpath) == "fast5":
#             for file in os.listdir(dirpath):
#                 if os.path.splitext(file)[1] == ".fast5":
#                     with open("guppy_inputs.txt", "a") as file:
#                         file.write(dirpath)
#                         file.write("\n")
#                     print(dirpath, "contains fast5 files")
#                     break
#             else:
#                 print(dirpath, "contains NO fast5 files")
