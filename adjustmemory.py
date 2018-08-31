#
# adjust the <job.memory> variable to rformat it for java -Xmx and reduce it by a couple of GB
# since if its set to "64 GB", then the docker container will get 64 GB but the java can't grab
# all that so needs to grab 2 GB less to leave some room for the OS  
#

import sys
val = sys.argv[1]
val = val.replace(" ","")
valInt = int(val) -2
val = str(valInt)

# unit could be G or Gb, strip the b if present
unit=sys.argv[2]
if len(unit) > 1:
    unit = unit[:-1]

print(val+unit)

