from ftplib import FTP_TLS
import logging
import threading
import socket
import os

ftpAddr ='cdaweb.gsfc.nasa.gov' #Address of the ftp
remoteDataDir = '/pub/data' #ftp database directory
verbose = True
downloadBlocksize = 1024*32

#loggingHandlers = [logging.FileHandler(os.path.join(logFileDir, 'download.log'))]
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
loggingHandlers = []
if verbose:
    loggingHandlers.append(logging.StreamHandler())
loggingHandlers[0].setFormatter(formatter)
logging.basicConfig(level=logging.DEBUG, handlers=loggingHandlers)

ftp = FTP_TLS(ftpAddr)
ftp.login()
ftp.sock.close()
ftp.getwelcome()
ftp.pwd()
ret = FTP_TLS().login(*argt, **argd)
