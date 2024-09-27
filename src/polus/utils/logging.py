import sys
import logging


global logger
logger = logging.getLogger(__name__)

def RaiseError(message):
    if (isinstance(message,str)):
        logger.error(msg=message)
        sys.exit()
    else:
        InvalidLogMessage()

def RaiseWarning(message):
    if (isinstance(message,str)):
        logger.warning(msg=message)
    else:
        InvalidLogMessage()

def PrintInfo(message):
    if (isinstance(message,str)):
        logger.info(msg=message)
    else:
        InvalidLogMessage()

def InvalidLogMessage():
    logger.warning(msg="Invalid logging message")
