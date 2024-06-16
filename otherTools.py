__author__ = 'Yufei Zhou'

## This file contain basic functions in python that require only low level models

from datetime import datetime
import copy
import numpy as np

##
class Directory:
    def __init__(self, dic=None, lis=None):
        '''
        Parameters:
            dic: or a dict, such as {'mms': {'mms1': {'fpi': {}}, 'mms2': {'fpi': {}, 'fgm': {}}}}
            lis: list of list of string, such as [['mms', 'mms1', 'fpi'], ['mms', 'mms2', 'fpi']]
        '''
        def addInfoToDic(dic):
            if isinstance(dic, dict):
                dic.update({'__info': {'type': 'directory'}})
                for key_, item_ in dic.items():
                    if key_ != '__info':
                        addInfoToDic(item_)
        self.dic = dic
        self.lis = lis

    def generate_dic_from_lis(self):
        '''
        Purpose:
            This function is to convert a list of list of str to a dict. Each leaf of the tree of the dictionary corresponds to an element in the list, in the form of [[key, key, key, ..., key, key], ...]
        '''
        assert self.dic is None and self.lis
        dic = {}
        ls = self.lis
        for l in ls:
            currentDic = dic
            for key in l:
                item = currentDic.get(key)
                if item is not None:
                    if isinstance(item, str):
                        break
                    elif isinstance(item, dict):
                        previousDic = currentDic
                        currentDic = item
                else:
                    currentDic[key] = {}
                    previousDic = currentDic
                    currentDic = previousDic[key]
        self.dic = dic

    def generate_lis_from_dic(self):
        assert self.lis is None and self.dic
        lis = dict2list(self.dic)
        self.lis = [ll[:-1] for ll in lis]

##
def dict2list(dic, withLeafDict=False):
    '''
    Purpose:
        This function is to convert a dictionary to a list. Each leaf of the tree of the dictionary corresponds to an element in the list
    Parameters:
        dic: a dictionary
        withLeafDict: if false the returned list is in the form of [[key, key, key, ..., key, value], ...]. Otherwise, the returned list is in the form of [[key, ..., key, key, value, dic[key][...][key]], ...]
    '''
    def f(dic, listTilNow_, finalList, withLeafDict):
        for key, value in dic.items():
            listTilNow = listTilNow_.copy()
            listTilNow.append(key)
            if isinstance(value, dict) and value:
                f(value, listTilNow.copy(), finalList, withLeafDict=withLeafDict)
            else:
                listTilNow.append(value)
                if withLeafDict:
                    listTilNow.append(dic)
                finalList.append(listTilNow)
    finalList = []
    listTilNow = []
    f(dic, listTilNow.copy(), finalList, withLeafDict)
    return finalList


def doAtLeavesOfADict(dic, do):
    '''
    Purpose:
        This function is to do something at each leaf of a dictionary
    Parameters:
        dic: a dictionary
        do: a function which has a list in the form of [key, key, key, ..., key, value, dic] as its parameter, where dic is the leaf dictionary
    '''
    lis = dict2list(dic, withLeafDict=True)
    rets = []
    for ll in lis:
        ret = do(ll)
        rets.append(ret)
    return rets

class DictTree(dict):
    def __init__(self, *args):
        super().__init__(*args)

    def toList(self, withLeafDict=False):
        finalList = dict2list(self, withLeafDict=withLeafDict)
        return finalList

    def atEachLeafDo(self, do):
        rets = doAtLeavesOfADict(self, do)
        return rets

    def difference(self, d2):
        selfSole = compareDictRecursively(self, d2)
        return DictTree(selfSole)

    def union(self, d2):
        u = unionDictRecursively(self, d2)
        return DictTree(u)

    def intersection(self, d2):
        its = intersectDictRecursively(self, d2)
        return DictTree(its)

    def print(self, depth=0):
        for key, item in self.items():
            print("\t-"*depth + key)
            if isinstance(item, dict):
                subDictTree = DictTree(item)
                subDictTree.print(depth=depth+1)

    def getSubDictTreeByKeys(self, keys):
        '''
        keys: a tuple or a list such as ['a', 'b', 'c']
        return:
            a dictTree obtained by using keys consecutively. For example, if keys = ['a', 'b', 'c'], return self['a']['b']['c']
        '''
        dic = self
        for key in keys:
            dic = dic[key]
        return DictTree(dic)

def compareDictRecursively(d1, d2):
    d1Sole = {}
    for key, d1Value in d1.items():
        try: d2Value = d2[key]
        except:
            d1Sole[key] = d1Value
            continue
        if d2Value == d1Value:
            continue
        elif isinstance(d2Value, dict):
            d1SoleKey = compareDictRecursively(d1Value, d2Value)
        else:
            d1SoleKey = d1Value
        d1Sole[key] = d1SoleKey
    return d1Sole


def intersectDictRecursively(d1, d2):
    intersection = {}
    for key, d1Value in d1.items():
        try: d2Value = d2[key]
        except: continue
        if d1Value == d2Value:
            its_ = d1Value
        elif d1Value is None:
            its_ = d2Value
        elif isinstance(d1Value, dict):
            if isinstance(d2Value, dict):
                its_ = intersectDictRecursively(d1Value, d2Value)
            elif d2Value is None:
                its_ = d1Value
            elif d2Value in d1Value.keys():
                its_ = {d2Value: d1Value[d2Value]}
            else:
                raise Exception('input error')
        elif d1Value in d2Value.keys():
            its_ = {d1Value: d2Value[d1Value]}
        else:
            raise Exception('input error')
        intersection[key] = its_
    return intersection


def unionDictRecursively(d1, d2, toRoot=True, priority=1):
    '''
    Parameters:
        toRoot: if true, the union is to the root of the two dict, e.g. d1={1:{2:None}} d2={1:{2:{3:None}}} gives {1:{2:None}}. If false it is to the leaves of the two dict, e.g. d1={1:{2:None}} d2={1:{2:{3:None}}} gives {1:{2:{3:None}}}
        priority: When d1 conflict with d2, this parameter determine how the handle it. Its possible values are 0, 1, and 2. Setting 0 to raise exception, 1 to choose d1, 2 to choose d2.
    '''
    union = copy.deepcopy(d2)
    for key, d1Value in d1.items():
        if key in d2:
            d2Value = d2[key]
            if d1Value is None:
                if toRoot:
                    union_ = d1Value
                else:
                    union_ = d2Value
            elif d1Value == d2Value:
                union_ = d1Value
            elif isinstance(d1Value, dict):
                if isinstance(d2Value, dict):
                    union_ = unionDictRecursively(d1Value, d2Value)
                elif d2Value is None:
                    if toRoot:
                        union_ = d2Value
                    else:
                        union_ = d1Value
                elif d2Value in d1Value.keys():
                    union_ = d1Value.copy()
                    if toRoot:
                        union_[d2Value] = None
                else:
                    raise Exception('input error')
            elif isinstance(d2Value, dict) and d1Value in d2Value.keys():
                union_ = d2Value.copy()
                if toRoot:
                    union_[d1Value] = None
            else:
                if priority == 0:
                    raise Exception('input error')
                elif priority == 1:
                    union_ = d1Value
                elif priority == 2:
                    union_ = d2Value
        else:
            union_ = d1Value
        union[key] = union_
    return union


def list2dict(ls):
    '''
    Purpose:
        This function is to convert a list of list of str to a dict. Each leaf of the tree of the dictionary corresponds to an element in the list, in the form of [[key, key, key, ..., key, value], ...]
    '''
    dic = {}
    for l in ls:
        currentDic = dic
        for key in l[:-1]:
            item = currentDic.get(key)
            if item is not None:
                if isinstance(item, str):
                    break
                elif isinstance(item, dict):
                    previousDic = currentDic
                    currentDic = item
            else:
                currentDic[key] = {}
                previousDic = currentDic
                currentDic = previousDic[key]
        else:
            previousDic[key] = l[-1]
    return dic


def datetime2list(dateTime, epochType='CDF_EPOCH'):
    if epochType == 'CDF_EPOCH':
        return [dateTime.year, dateTime.month, dateTime.day, dateTime.hour, dateTime.minute, dateTime.second, dateTime.microsecond//10**3]
    if epochType == 'CDF_EPOCH16':
        return [dateTime.year, dateTime.month, dateTime.day, dateTime.hour, dateTime.minute, dateTime.second, dateTime.microsecond//10**3, dateTime.microsecond%1000, 0, 0]
    elif epochType == 'CDF_TIME_TT2000':
        return [dateTime.year, dateTime.month, dateTime.day, dateTime.hour, dateTime.minute, dateTime.second, dateTime.microsecond//10**3, dateTime.microsecond%1000, 0]


def decimal2binaryArray(array, order='>'):
    binary = []
    while np.any(array > 0):
        array, rem = np.divmod(array, 2)
        binary.append(rem)
    if order == '>':
        binary.reverse()
    elif order == '<':
        pass
    binary_array = np.array(binary)
    return np.moveaxis(binary_array, 0, -1)


def sizeof_fmt(num, suffix="B", fmt=True):
    if fmt:
        for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
            if abs(num) < 1024.0:
                return f"{num:3.1f}{unit}{suffix}"
            num /= 1024.0
        return f"{num:.1f}Yi{suffix}"
    else:
        for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
            if abs(num) < 1024.0:
                return f"{num}{unit}{suffix}"
            num /= 1024.0
        return f"{num}Yi{suffix}"

def sizeof(string):
    '''
    transform 2.3M into 2.3*1024**2

    '''
    unit_map = {'K': 1, 'M': 2, 'G': 3, 'T': 4, 'P': 5, 'E': 6, 'Z': 7}
    try:
        _ = int(string[-1])
        return float(string)
    except:
        return float(string[:-1]) * 1024**unit_map[string[-1]]

def round_number(num, ndigs=0):
    if ndigs:
        return round_number(num*10**ndigs)/10**ndigs
    if num - int(num) == 0.5:
        return int(np.ceil(num))
    else:
        return round(num)
