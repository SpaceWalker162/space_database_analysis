import databaseUserTools as dut
from datetime import datetime

databaseDir = '/home/yufei/Documents/remoteDatabase/data'
workDataDir = '/home/yufei/Documents/remoteDatabase/data'
mission = 'cluster'
spacecraftNames = ['c1', 'c2', 'c3', 'c4']
datasets = {mission: {}}
for spacecraftName in spacecraftNames:
    instrumentation = spacecraftName.upper() + '_CP_FGM_5VPS'
    datasets[mission].update({spacecraftName: instrumentation})
print(datasets)
interval = (datetime(2001, 1, 1), datetime(2020, 1, 1))
dut.extractFiles(databaseDir, workDataDir, datasets, interval, keepIfExist=False)
