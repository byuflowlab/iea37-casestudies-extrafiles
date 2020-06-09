#### TO DO ####
# [ ] Get read in of boundary written
# [ ] Get splining of boundary with flowmath written
# [ ] Get checkBoundary written

import YAML
using Printf
using Parameters

function getBoundaryCs4YAML(file_name)
    ### Retrieve boundary inflection points for IEA37 case study 4 from the <.yaml> file.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    bndrs = f["boundaries"]

    ptList3a = bndrs["IIIa"]
    ptList3b = bndrs["IIIb"]
    ptList4a = bndrs["IVa"]
    ptList4b = bndrs["IVb"]
    ptList4c = bndrs["IVc"]

    return ptList3a, ptList3b, ptList4a, ptList4b, ptList4c
end

fn = "../startup-files/iea37-boundary-cs4.yaml"

#ptsBndry = getBoundaryCs4YAML(fn)
plL3a, plL3b, plL4a, plL4b, plL4c = getBoundaryCs4YAML(fn)
# plL3a = ptsBndry[1]
# plL3b = ptsBndry[2]
# plL4a = ptsBndry[3]
# plL4b = ptsBndry[4]
# plL4c = ptsBndry[5]

print(typeof(plL3a))