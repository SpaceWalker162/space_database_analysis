import os
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
        cdf_master = cdflib.CDF(src)
        if (cdf_master.file != None):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            info = cdf_master.cdf_info() # Get the cdf's specification
            cdf_file = cdflib.cdfwrite.CDF(dst, delete=True)
            globalAttrs = cdf_master.globalattsget() # Get the global attributes
            # For the moment there is a bug in the cdflib.cdfread.CDF.globalattsget(). The following lines tackle this.
            for att, item in globalAttrs.items():
                dic = {}
                for i, it in enumerate(item):
                    dic[i] = it
                globalAttrs[att] = dic
            cdf_file.write_globalattrs(globalAttrs)
            zvars = info.zVariables
            for x in range (0, len(zvars)):
                # For the moment there is this bug in varinq(), which returns a VDRInfo object rather than a dictionary to be accepted by write_var(). Therefore we use __dict__
                varattrs = cdf_master.varattsget(zvars[x])
                varinfo = cdf_master.varinq(zvars[x]).__dict__
                try:
                    vardata = cdf_master.varget(zvars[x])
                except:
                    vardata = None
                varname = varinfo['Variable']
                if varinfo['Data_Type'] not in [31, 32, 33]:
                    if varname[-4:] == 'time':
                        varname = varname[:-4] + 'epoch'
                        dim_size = vardata.shape[1:]
                        varinfo = {
                                'Variable': varname,
                                'Data_Type': 33,
                                'Num_Elements': 1,
                                'Rec_Vary': True,
                                'Dim_Sizes': dim_size,
                                'Block_Factor': 5462,
                                    }
                        varattrs = {
                                'CATDESC': varname + ' (CDF_TIME_TT2000)',
                                'FIELDNAM': varname,
                                'FILLVAL': -1e+31,
                                'VALIDMIN': 59958200000000.0,
                                'VALIDMAX': 66301200000000.0,
                                'VAR_TYPE': 'support_data',
                                'LABLAXIS': 'UT',
                                }
                        vardata = cdflib.epochs.CDFepoch.timestamp_to_tt2000(vardata)
                    else:
                        if 'DEPEND_0' in varinfo:
                            varinfo['DEPEND_0'] = '_'.join(varinfo['DEPEND_TIME'].split('_')[:-1] + ['epoch'])
                            del varinfo['DEPEND_TIME']
                            del varinfo['DEPEND_EPOCH0']
                else:
                    assert vardata is None or varname == 'range_epoch' or vardata == 62167219200000.0
                    continue

                if varinfo['Data_Type'] in [31, 32, 33]:
                    varinfo['Compress'] = 0
                else:
                    varinfo['Compress'] = 6

                cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
            cdf_file.close()

