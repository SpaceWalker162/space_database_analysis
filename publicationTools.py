import re
import os
import subprocess


def extract_environement_from_tex(string, environmentName):
    texEnvironmentPattern = '\\\\begin{{{environmentName}}}.*?\\\end{{{environmentName}}}'.format(environmentName=environmentName)
    texEnvironmentPatternObj = re.compile(texEnvironmentPattern, flags=re.DOTALL)
    environmentStringList = texEnvironmentPatternObj.findall(string)
    return environmentStringList


def findMatchingDelimiter(string, start=0, left='{', right='}'):
    '''
    Example: string='abc{efl{gl}}', start=4, return 12
    '''
    leftCount = 1
    for sInd, s in enumerate(string[start:]):
        if s == left:
            leftCount += 1
        elif s == right:
            leftCount -= 1
        if leftCount == 0:
            return sInd + start

def extract_caption_from_environment(environmentString):
    '''
    Purpose: extract caption from environment string in tex
    '''
    captionNameString = r'\caption{'
    sLen = len(captionNameString)
    startInd = environmentString.find(captionNameString) + sLen
    endInd = findMatchingDelimiter(string=environmentString, start=startInd)
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
    stringInHTML4Format = pandoc_transform(string, oldFormat='latex', newFormat='html4')
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
    finalString = ''.join(processedStringList)
    return finalString
