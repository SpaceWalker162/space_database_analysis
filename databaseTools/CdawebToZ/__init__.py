import os
import numpy as np
import cdflib


class CdawebTHEMISFile:

    def __init__(self, path=None):
        self.path = path

    @staticmethod
    def convert_from_cdaweb_to_z(src, dst):
        '''
        epoch data compress level is 0 by default. Other data's compress level is 6 by default.
        Parameters:
            src: str. Source file path.
            dst: str. Destination file path.
        '''
        try:
            cdf_master = cdflib.CDF(src)
        except:
            return
        if (cdf_master.file != None):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            info = cdf_master.cdf_info() # Get the cdf's specification
            globalAttrs = cdf_master.globalattsget() # Get the global attributes
            # For the moment there is a bug in the cdflib.cdfread.CDF.globalattsget(). The following lines tackle this.
            for att, item in globalAttrs.items():
                dic = {}
                for i, it in enumerate(item):
                    dic[i] = it
                globalAttrs[att] = dic
            zvars = info.zVariables
            vardata_list = []
            varattrs_list = []
            varinfo_list = []
            for x in range (0, len(zvars)):
                # For the moment there is this bug in varinq(), which returns a VDRInfo object rather than a dictionary to be accepted by write_var(). Therefore we use __dict__
                varattrs = cdf_master.varattsget(zvars[x])
                varinfo = cdf_master.varinq(zvars[x]).__dict__
                if varinfo['Last_Rec'] >= 0:
                    vardata = cdf_master.varget(zvars[x])
                else:
                    vardata = None
                varname = varinfo['Variable']
                if varinfo['Data_Type'] not in [31, 32, 33]:
                    if varname[-4:] == 'time':
                        varname = varname[:-4] + 'epoch'
                        varinfo['Variable'] = varname
                        varinfo['Data_Type'] = 33
#                        varinfo = {
#                                'Variable': varname,
#                                'Data_Type': 33,
#                                'Num_Elements': 1,
#                                'Rec_Vary': True,
#                                'Dim_Sizes': dim_size,
#                                'Block_Factor': 5462,
#                                'Sparse': 'No_sparse'
#                                    }
                        varattrs = {
                                'CATDESC': varname + ' (CDF_TIME_TT2000)',
                                'FIELDNAM': varname,
                                'FILLVAL': -1e+31,
                                'VALIDMIN': 59958200000000.0,
                                'VALIDMAX': 66301200000000.0,
                                'VAR_TYPE': 'support_data',
                                'LABLAXIS': 'UT',
                                }
                        if varinfo['Last_Rec'] >= 0:
                            vardata = cdflib.epochs.CDFepoch.timestamp_to_tt2000(vardata)
                    else:
                        if 'DEPEND_0' in varattrs:
                            varattrs['DEPEND_0'] = '_'.join(varattrs['DEPEND_TIME'].split('_')[:-1] + ['epoch'])
                            del varattrs['DEPEND_TIME']
                            del varattrs['DEPEND_EPOCH0']
                else:
                    assert vardata is None or varname == 'range_epoch' or vardata == 62167219200000.0
                    continue

                if varinfo['Data_Type'] in [31, 32, 33]:
                    varinfo['Compress'] = 0
                else:
                    varinfo['Compress'] = 6
                varinfo_list.append(varinfo)
                varattrs_list.append(varattrs)
                vardata_list.append(vardata)
            cdf_file = cdflib.cdfwrite.CDF(dst, delete=False)
            cdf_file.write_globalattrs(globalAttrs)
            for varinfo, varattrs, vardata  in zip(varinfo_list, varattrs_list, vardata_list):
                if varinfo['Sparse'].lower() == 'no_sparse':
                    cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
                elif varinfo['Sparse'].lower() == 'prev_sparse':
                    prevrecs = np.insert((np.nonzero(np.diff(vardata, axis=0))[0] + 1), 0, 0)
                    cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=[prevrecs, vardata])
            cdf_file.close()

