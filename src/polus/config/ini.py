#!/usr/bin/env python
import os
import sys
from configparser import ConfigParser

global keys
keys                = ["atoms","name","props","natoms","method","train","val","test","iqa_filt","q00_filt"]
mandatory_sections  = ["system","sampling","filtering"]

def check_config(root: str) -> None:
    config_file = os.path.join(root,"polus.ini")
    if os.path.isfile(config_file):
        config=ConfigParser()
        config.read(config_file)
        config_sections = set(config.sections())
        if (config_sections!=set(mandatory_sections)):
            print("POLUS| Program complains::: Missing sections in polus.ini file")
            sys.exit()
        else:
            keys1 = config["system"]
            keys2 = config["sampling"]
            keys3 = config["filtering"]
            config_keys = list(set(keys1) | set(keys2) | set(keys3))
            for key in keys:
                if (key not in config_keys):
                    print("POLUS| Program complains::: Missing key [{}] in polus.ini file".format(key))
                    sys.exit()
    else:
        print("POLUS| Program complains::: {} file not found".format(config_file))
        sys.exit()
