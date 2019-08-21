#
# adjust the <job.memory> variable to rformat it for java -Xmx and reduce it by a couple of GB
# since if its set to "64 GB", then the docker container will get 64 GB but the java can't grab
# all that so needs to grab 1 GB less to leave some room for the OS  
#

import sys
val = sys.argv[1]
val = val.replace(" ","")
try:
    valInt = int(val) -1
    val = str(valInt)
    # unit could be G or Gb, strip the b if present
    unit=sys.argv[2]
    if len(unit) > 1:
        unit = unit[:-1]

    print(val+unit)

except:
    # probably came in as '16G' or '16Gb' without a space
    # just need to strip the 'b' if it ends with one
    if (val[-1].lower() == 'b'):
        print(val[:-1])
    else:
        print(val)


