# Handle multiple (how? -- split them maybe?)
#  Find the approprite file for the run to use for QC (seq summ?)
import os
def FindFast5(run_path_list):
    fast5_path_list=[]
    count = 0
    for run_path in run_path_list:
        for dirpath, dirnames, filenames in os.walk(run_path, topdown = False):
            print(dirpath, count)
            count = count + 1
            if os.path.basename(dirpath) == "fast5":
                for file in os.listdir(dirpath):
                    if os.path.splitext(file)[1] == ".fast5":
                        fast5_path_list.append(dirpath)
                    else:
                        print(dirpath, "contains NO fast5 files")
    return fast5_path_list



FindFast5(['test/temp1', 'test/temp2', 'test/temp3', 'test/temp4'])
