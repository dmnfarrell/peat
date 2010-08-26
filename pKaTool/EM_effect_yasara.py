# YASARA PLUGIN
# TOPIC:       pKaTool/EM_effect titration interface
# TITLE:       Start pKaTool or EM_effect
# AUTHOR:      Jens Nielsen
# LICENSE:     GPL
# DESCRIPTION: This plugin starts EM_effect or pKaTool
#
# This is a YASARA plugin to be placed in the /plg subdirectory
# Go to www.yasara.com/plugins for documentation and downloads
#
"""
MainMenu after Analyse: _p_Ka/NMR
  PullDownMenu : _p_KaTool
    Request: StartpKaTool
  PullDownMenu : _E_Meffect
    Request: StartEMeffect
"""
#
# Add the Yasara directory to the search path
#
import sys,os
script_location=sys.argv[1]
plgdir=os.path.split(script_location)[0]
yasara_dir=os.path.split(plgdir)[0]
sys.path.append(yasara_dir)
sys.path.append(plgdir)
#
# Import Yasara
#
import yasara,string,ftplib,urllib,disk,time
#
# Get the location of the pKaTool/EM_effect
#
pKadir='/home/people/nielsen/lib/development/pKaTool/'
sys.path=[pKadir]+sys.path

if (yasara.request=="StartpKaTool"):
  import pkaTool
  pKaTool
elif (yasara.request=="StartEMeffect"):
  import EM_effect
  EM_effect.EM_effect(Yasara=yasara).mainloop()

# THIS MUST ALWAYS BE THE LAST COMMAND
yasara.plugin.end()
