import databaseUserTools as dut
from datetime import datetime

databaseDir = '/media/yufei/server'
workDataDir = '/media/yufei/server'
mission = 'Cluster'
spacecraftNames = ['C1', 'C2', 'C3', 'C4']
datasets = {mission: {}}
for spacecraftName in spacecraftNames:
    instrumentation = spacecraftName + '_CP_FGM_FULL'
    datasets[mission].update({spacecraftName: instrumentation})
print(datasets)
interval = (datetime(2001, 1, 1), datetime(2010, 1, 1))
dut.extractFiles(databaseDir, workDataDir, datasets, interval, keepIfExist=False)
