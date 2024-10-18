import sys
import logging


global logger
logger = logging.getLogger(__name__)

def RaiseError(message: str) -> None:
    if (isinstance(message,str)):
        logger.error(msg=message)
        sys.exit()
    else:
        InvalidLogMessage()

def RaiseWarning(message: str) -> None:
    if (isinstance(message,str)):
        logger.warning(msg=message)
    else:
        InvalidLogMessage()

def PrintInfo(message: str) -> None:
    if (isinstance(message,str)):
        logger.info(msg=message)
    else:
        InvalidLogMessage()

def InvalidLogMessage() -> None:
    logger.warning(msg="Invalid logging message")
