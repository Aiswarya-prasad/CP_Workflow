import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

counts = {}
total = 0
with open(snakemake.input[0], 'r') as tsv, open(snakemake.output[0], 'w') as out:
    header = tsv.readline().split()
    for line in tsv.readlines():
        s = line.strip().split()
        total += 1
        barcode = s[2]
        if barcode == "none":
            pass
        else:
            if int(barcode) < 10:
                barcode = "0" + barcode
            # barcode = "barcode" + barcode
        if barcode in counts.keys():
            counts[barcode] += 1
        else:
            counts.__setitem__(barcode, 0)
    lines = []
    for barcode_num, count in counts.items():
        lines.append("{}: {} ({}%)\n".format(barcode_num, count, round(count/total*100, 2)))
    out.writelines(lines)

r = range(len(counts.keys()))
rawData = {"green": [], "orange": []}
percentages = {}
for item in sorted(counts):
    rawData["green"].append(counts[item])
    rawData["orange"].append(total - counts[item])
df = pd.DataFrame(rawData)
totals = [i+j for i,j in zip(df["green"], df["orange"])]
green = [i / j * 100 for i,j in zip(df['green'], totals)]
orange = [i / j * 100 for i,j in zip(df['orange'], totals)]

# plot
barH = 1
names = sorted(counts)
# Create green Bars
# Create orange Bars
rect = plt.barh(r, green, color='#b5ffb9', edgecolor='white', height=barH)
plt.barh(r, orange, left=green, color='#f9bc86', edgecolor='white', height=barH)
plt.yticks(r, names)
plt.ylabel("Barcode")
barNum = 0
for box in rect:
    width = box.get_width()
    if width > 1:
        plt.text(width + 2, box.get_y() + box.get_height()/2, "{}".format(counts[sorted(counts)[barNum]]), ha='left', va='center')
    barNum += 1
plt.savefig(snakemake.ouput[1])
