import sys
import os
import space_database_analysis.databaseTools.CdawebToZ as CTZ
import space_database_analysis.databaseTools as dbt

if __name__ == '__main__':
    # srcdata and dstdata are the paths to the source/destination data directories. For example dstdata might be /mnt/pub/data
    srcdata, dstdata = sys.argv[1:3]
    srcdata = os.path.expanduser(srcdata)
    dstdata = os.path.expanduser(dstdata)

    verbose = True
    srcdataTree = dbt.readFileInfoRecursively(path=srcdata, verbose=verbose, facts='stats')
    toD = srcdataTree.toList()
    toDNew = []
    for toD_ in toD:
        if 'fgm' in toD_ and 'l2' in toD_:
            toDNew.append(toD_)
    toD = toDNew
    total_counts = len(toD)
    for ind, src_name_in_list in enumerate(toD):
        src = os.path.join(srcdata, *src_name_in_list[:-1])
        ss_, ext = os.path.splitext(src)
        if ext.lower() == '.cdf':
            ss_, ext = os.path.splitext(ss_)
            if ext == '.z':
                continue
            else:
                pass
        else:
            continue
        dst = os.path.splitext(os.path.join(dstdata, *src_name_in_list[:-1]))[0] + '.z.cdf'
        if os.path.exists(dst):
            continue
#        if src_name_in_list[0] == 'themis' and 'fgm' in src_name_in_list:
        if verbose:
            print('converting ', src)
        CTZ.CdawebTHEMISFile.convert_from_cdaweb_to_z(src, dst)
        print('done, progress {}/{}'.format(ind, total_counts))

#import cdflib
#filename = os.path.expanduser('~/Documents/remoteDatabase/data/themis/tha/l2/fgm/2008/tha_l2_fgm_20080109_v01.cdf')
#cdfFile = cdflib.CDF(filename)
#cdfInfo = cdfFile.cdf_info()
#
#cdfFile.varinq('tha_fgs_epoch')
#cdfFile.varget('tha_fgs_epoch0')
#cdfFile.varinq('tha_fgs_time')
