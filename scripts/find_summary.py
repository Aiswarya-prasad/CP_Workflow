# Handle multiple (how? -- split them maybe?)
#  Find the approprite file for the run to use for QC (seq summ?)
#
# EDIT the function below to find required files any time
# disabling topdown starts with directories at the bottom (nodes with single leaf)
# 
# def FindFast5(run_path_list):
#     """
#     when provided with a list of directory paths with
#     raw data (per run), it returns a list of paths to
#     the fast5 directories in them which contain fast5 files
#
#     it can handle multiple fast5 directories in each
#     parent given
#     the path will contain the parent directory (run name)
#     """
#     # for root, dirs, files in os.walk(run_path, topdown=False):
#     #     for name in dirs:
#     #         if os.path.basename(dirpath) == "fast5":
#     #             print(name)
#
#     fast5_path_list=[]
#     for run_path in run_path_list:
#         for dirpath, dirnames, filenames in os.walk(run_path, topdown = False):
#             if os.path.basename(dirpath) == "fast5":
#                 for file in os.listdir(dirpath):
#                     if os.path.splitext(file)[1] == ".fast5":
#                         fast5_path_list.append(dirpath)
#                     else:
#                         print(dirpath, "contains NO fast5 files")
#     return fast5_path_list
