#-*-Python-*-
# Sconstruct for nut
# T. M. Kelley
# Dec 18, 2010
# Input to scons
# (c) Copyright 2010 LANSLLC all rights reserved.

import os

# external pieces:
scons_dir = "./scons"
execfile(scons_dir + "/scons_master")
execfile(scons_dir + "/scons_project")
execfile(scons_dir + "/random123")

# user can influence preprocessing from environment 
try:
    external_cppflags = os.environ['CPPFLAGS']
    cppflags += " %s " % external_cppflags
except KeyError: pass

# load the build environment
compiler_file = "scons_compiler_%s_%s_%s" % (par, arch, comp_family)
execfile(scons_dir + "/%s" % compiler_file)

Export('env','build_root','export_libdir','export_incdir','export_bindir',
       'project_incdir','project_libdir',"libs","libpaths",
       "cpppaths","cxxflags","cflags","cppflags","cppdefines",
       "ldflags","arflags"
       )

# recurse into lib, test, and app directories
SConscript(
    ["lib/SConscript",
     "apps/SConscript",
     "test/SConscript",
     ])

# version
__id__ = "$Id$"

# End of file

