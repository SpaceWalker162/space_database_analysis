import os
import numpy as np
import cdflib
from cdflib import cdfwrite, cdfread
import subprocess as sp

def convert_cluster_file(src_file, dst_file):
    cdf_master = cdfread.CDF(src_file)
    if (cdf_master.file != None):
        # Get the cdf's specification
        info = cdf_master.cdf_info()
        cdf_file = cdfwrite.CDF(dst_file, delete=True)

        # Get the global attributes
        globalAttrs = cdf_master.globalattsget()
        for att, item in globalAttrs.items():
            dic = {}
            for i, it in enumerate(item):
                dic[i] = it
            globalAttrs[att] = dic
        # Write the global attributes
        cdf_file.write_globalattrs(globalAttrs)
        zvars = info.zVariables
    #    print('no of zvars=',len(zvars))
        # Loop thru all the zVariables
        for x in range (0, len(zvars)):
            # Get the variable's specification
            varinfo = cdf_master.varinq(zvars[x]).__dict__
            varinfo['Compress'] = 6
            #print('Z =============>',x,': ', varinfo['Variable'])
            # Get the variable's attributes
            varattrs = cdf_master.varattsget(zvars[x])
            if (varinfo['Sparse'].lower() == 'no_sparse'):
                # A variable with no sparse records... get the variable data
                vardata = cdf_master.varget(zvars[x])
                # Create the zVariable, write out the attributes and data
                cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
            else:
                # A variable with sparse records...
                vardata = cdf_master.varget(zvars[x])
                # Create the zVariable, write out the attributes and data
                cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
        cdf_file.close()

def convert_under_dir(path_to_dir, path_to_dst_dir):
    for file in os.scandir(path_to_month):
        if file.is_file():
            src_file = os.path.join(path_to_dir, file.name)
            datasetID, timeEXT = file.name.split('__')
            timeEXT_list = timeEXT.split('_')
            tStr = timeEXT_list[0]
            versionStr = 'v20' + timeEXT_list[-1][1:]
            dst_file = os.path.join(path_to_dst_dir, '_'.join([datasetID.lower(), tStr, versionStr]))
            convert_cluster_file(src_file, dst_file)

clusterPath = '/mnt/pub/data/cluster'

clusterInds = range(1, 5)
for i in range(1, 5):
    clusterName = 'c' + str(i)
    clusterIndPath = os.path.join(clusterPath, clusterName)
    os.chdir(clusterIndPath)
    path_from_tos = [('C{}_CP_FGM_FULL'.format(i), 'fgm/mfield_3dvect_fullreso'), ('C{}_CP_FGM_5VPS'.format(i), 'fgm/mfield_3dvect_5vectpersec')]
    for path_from, path_to in path_from_tos:
        for yearD in os.scandir(path_from):
            if yearD.is_dir():
                path_to_year = os.path.join(path_from, yearD.name)
                print('in {} {}'.format(clusterName, path_to_year))
                path_to_dst_dir = os.path.join(path_to, yearD.name)
                os.makedirs(path_to_dst_dir, exist_ok=True)
                if 'FULL' in path_from:
                    for monthD in os.scandir(path_to_year):
                        if monthD.is_dir():
                            path_to_month = os.path.join(path_to_year, monthD.name)
                            convert_under_dir(path_to_month, path_to_dst_dir)
                elif '5VPS' in path_from:
                    convert_under_dir(path_to_year, path_to_dst_dir)
