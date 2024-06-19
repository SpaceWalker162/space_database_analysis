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
    for src_name_in_list in toD:
        src = os.path.join(srcdata, *src_name_in_list[:-1])
        dst = os.path.splitext(os.path.join(dstdata, *src_name_in_list[:-1]))[0] + '.z.cdf'
#        if src_name_in_list[0] == 'themis' and 'fgm' in src_name_in_list:
        if 'fgm' in src_name_in_list and 'l2' in src_name_in_list:
            if verbose:
                print('converting ', src)
            CTZ.CdawebTHEMISFile.convert_from_cdaweb_to_z(src, dst)

