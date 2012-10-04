#!/usr/bin/env python
# gen_config.py
# T. M. Kelley
# Dec 18, 2010
# (c) Copyright 2010 LANSLLC all rights reserved.

import time,os


def generate_config_h( target, source, env):
    """generate_config_h( target, source, env)
    generate the config.h file for nut.
    Parameters are supplied by SCons.

    The target file will be config.h, the source file is config.h.pyin.
    Parameters are passed in from SConstruct/SConscript through the
    build environment's 'config' dictionary; these are used to replace
    text in the source file (config.h.pyin) to produce the text to put
    into the target file (config.h)
    """

    intext = file(source[0].path).read()
    outtext = intext

    config = env['config']

    ## Example use:
    # value = config['key']
    # outtext = outtext.replace("<some key>","#define SOMEKEY value")

    # substitute date & year
    date = time.strftime("%b %d, %Y")
    year = time.strftime("%Y")
    outtext = outtext.replace("<date>", date)
    outtext = outtext.replace("<year>", year)

    outfile = file(target[0].path,'w')
    outfile.write(outtext)
    outfile.close()
    return


# version
__id__ = "$Id: template_gen_config.py,v 1.1.1.1 2008/12/02 23:27:42 tkelley Exp $"

# End of file
