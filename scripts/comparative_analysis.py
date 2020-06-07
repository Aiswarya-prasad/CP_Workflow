ListOfExpectedBarcodes = []
confDict = sample_dict:
    'Run0':
      '04': '01'
    'Run1':
      '02': '03'
      '03': '13'
    'Run2':
      '04': '06'
      '05': '07'
      '06': '10'
      '07': '14'
    'Run3':
      '08': '11'
      '09': '12'
      '10': '15'
      '11': '17'
      '12': '18'
    'Run4':
      '01': '19'
      '02': '20'
      '03': '21'
      '04': '22'
for RunName in confDict:
    for barCode in confDict[RunName].keys():
        ListOfExpectedBarcodes.append(join("qcat_trimmed", RunName, "barcode"+barCode+".fastq"))
