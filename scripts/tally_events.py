#!/wsu/el7/pre-compiled/python/3.7/bin/python
import glob
import sys

data = dict()
print(sys.argv)

which="HT"
if len(sys.argv)>1:
    which=sys.argv[1]


triggers = dict()
for file in glob.glob(f'out/trigger_count/{which}/pAu*.list'):
    for line in open(file).readlines()[1:]:
        vals = line.split()

        run = int(vals[0])
        trigger = int(vals[1])
        nevents = int(vals[2])

        if run not in data:
            data[run] = dict()

        if trigger in data[run]:
            data[run][trigger] += nevents
        else:
            data[run][trigger] = nevents

        if trigger in triggers:
            triggers[trigger] += nevents
        else:
            triggers[trigger] = nevents


with open(f"out/trigger_count/{which}/runid.list",'w') as fout:
    for runid in sorted(data.keys()):
        fout.write(f'{runid}\n')

with open(f"out/trigger_count/{which}/trigger.list",'w') as fout:
    for trigger in sorted(triggers.keys()):
        fout.write(f'{trigger:8} {triggers[trigger]:10}\n')

with open(f"out/trigger_count/{which}/trigger_cnt.list",'w') as fout:
    for runid in sorted(data.keys()):
        for trigger in sorted(data[runid].keys()):
            fout.write(f'{runid} {trigger:6} {data[runid][trigger]}\n')
