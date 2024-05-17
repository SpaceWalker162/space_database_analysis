import re
import os
import subprocess
import logging


class TexFile:
    def __init__(self, filePath=None, content=None):
        if filePath:
            self.filePath = filePath
            with open(self.filePath, 'r') as f:
                texContent = f.read()
            self.content = texContent
        elif content:
            self.content = content
        if self.content:
            documentStartInd = self.content.find(r'\begin{document}')
            self.preamble = self.content[:documentStartInd]
            self.document_content = self.content[documentStartInd:]

    def extract_environement(self, environmentName=None):
        return extract_environement_from_tex(string=self.document_content, environmentName=environmentName)

    def extract_macro(self, macroName=None):
        return extract_macro_from_tex(self.document_content, macroName=macroName)

def makeDelimiterPairDict(delimiterPairs):
    dic = {}
    for pair in delimiterPairs:
        if pair[0] == pair[1]:
            dic.update({pair[0]: {'s': pair[0], 'loc': None}})
            dic[pair[0]].update({'pair': dic[pair[0]]})
        else:
            dic.update({pair[0]: {'s': pair[0], 'loc': 'left'}})
            dic.update({pair[1]: {'s': pair[1], 'loc': 'right'}})
            dic[pair[0]].update({'pair': dic[pair[1]]})
            dic[pair[1]].update({'pair': dic[pair[0]]})
    return dic


def extract_environement_from_tex(string, environmentName):
    texEnvironmentPattern = '\\\\begin{{{environmentName}}}.*?\\\end{{{environmentName}}}'.format(environmentName=environmentName)
    texEnvironmentPatternObj = re.compile(texEnvironmentPattern, flags=re.DOTALL)
    environmentStringList = texEnvironmentPatternObj.findall(string)
    return environmentStringList

def extract_macro_from_tex(string, macroName=None):
    r'''
    Parameters:
        macroName: if None, extract all macros. If a string such as ref, extract \ref plus all braces following it.
    '''
    macroRanges = []
    macroContents = []
    if macroName is None:
        patStr = r'\[a-zA-Z]+'
        pat = re.compile(patStr)
    elif isinstance(macroName, str):
        patStr = '\\\\' + macroName
        pat = re.compile(patStr)
    matchCount = 0
    for match in pat.finditer(string):
        matchCount += 1
        start = match.end()
        immediateC = string[start]
        while immediateC in ['{', '<', '[']:
            right = delimiterPairDict[immediateC]['pair']['s']
            delimiterPairEnd = findMatchingDelimiter(string, start=start, left=immediateC, right=right)
            start = delimiterPairEnd+1
            immediateC = string[start]
        macroRange = (match.start(), delimiterPairEnd+1)
        macroRanges.append(macroRange)
        macroContent = string[slice(*macroRange)]
        macroContents.append(macroContent)
    return macroRanges, macroContents


def findMatchingDelimiter(string, start=0, left='{', right='}'):
    '''
    Example: string='abc{efl{gl}}', start=4, return 12
    '''
    logging.debug('fingMatchingDelimiter left: '+left)
    logging.debug('fingMatchingDelimiter right: '+right)
    if string[start] == right:
        initialLeftNetCount = -1
        stringToSearch = string[start::-1]
    elif string[start] == left:
        initialLeftNetCount = 1
        stringToSearch = string[start:]
    else:
        return None
    leftNetCount = 0
    for sInd, s in enumerate(stringToSearch):
        if s == left:
            leftNetCount += 1
        elif s == right:
            leftNetCount -= 1
        if leftNetCount == 0:
            return initialLeftNetCount*sInd + start

def extract_caption_from_environment(environmentString):
    '''
    Purpose: extract caption from environment string in tex
    '''
    captionNameString = r'\caption{'
    sLen = len(captionNameString)
    startInd = environmentString.find(captionNameString) + sLen
    endInd = findMatchingDelimiter(string=environmentString, start=startInd-1)
    caption = environmentString[startInd:endInd]
    return caption


def extract_figure_caption_from_tex(string):
    environmentStringList = extract_environement_from_tex(string, environmentName='figure')
    captions = []
    for environmentString in environmentStringList:
        caption = extract_caption_from_environment(environmentString)
        captions.append(caption)
    return captions

def pandoc_transform(string, oldFormat='latex', newFormat='html4'):
        commandString = 'pandoc --from {oldFormat} --to {newFormat} --no-highlight'.format(oldFormat=oldFormat, newFormat=newFormat)
        command = commandString.split(' ')
        process = subprocess.run(command, stdout=subprocess.PIPE, input=string, encoding='utf-8', capture_output=False)
        stringInNewFormat = process.stdout
        return stringInNewFormat


def nature_communications_figure_caption_string_transform(string):
    '''
    Parameters:
        repl: replace macro found in the string with some predefined content. For example if repl='[PLACEHOLDER#]', the first found macro will be replaced with [PLACEHOLDER1], the second with [PLACEHOLDER2], etc.
    '''
    macroRanges, macroContents = extract_macro_from_tex(string, macroName='cite')
    repl = '[PLACEHOLDER#]'
    addLen = 0
    for macroInd, macroRange in enumerate(macroRanges):
        macroRepl = repl.replace('#', str(macroInd+1))
        macroReplLen = len(macroRepl)
        macroLen = macroRange[1] - macroRange[0]
        start = macroRange[0] + addLen
        end = macroRange[1] + addLen
        string = string[:start] + macroRepl + string[end:]
        addLen += macroReplLen - macroLen
    stringInHTML4Format = pandoc_transform(string, oldFormat='latex', newFormat='html4')
#    logging.warning('\cite is lost in the conversion')
    stringInHTML4Format = stringInHTML4Format.replace('\n', ' ')

    htmlTagPattern = '(<.*?>)'
    htmlTagPatternObj = re.compile(htmlTagPattern)
    emPattern = '<(/?)em>'
    emPatternObj = re.compile(emPattern)
    subPattern = '<(/?)sub>'
    subPatternObj = re.compile(subPattern)

    stringList = htmlTagPatternObj.split(stringInHTML4Format)
    processedStringList = []
    for subInd, subString in enumerate(stringList):
        if subInd % 2 == 1:
            if match := emPatternObj.fullmatch(subString):
                processedString = "<{}i>".format(match.group(1))
                processedStringList.append(processedString)
            elif match := subPatternObj.fullmatch(subString):
                processedString = subString
                processedStringList.append(processedString)
            match = None
        else:
            processedStringList.append(subString)
    string = ''.join(processedStringList)

    for macroInd, macroContent in enumerate(macroContents):
        macroRepl = repl.replace('#', str(macroInd+1))
        string = string.replace(macroRepl, macroContent)
    return string

delimiterPairs = ['{}', '[]', '<>', '!!', '||', '$$', '&&']
delimiterPairDict = makeDelimiterPairDict(delimiterPairs)
