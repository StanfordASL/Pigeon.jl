import yaml
import glob
from osprey.msg import path
world2path = {
    "AxDes_mps2": "Ax_des_mps2",
    "edgeL_m": "edge_L_m",
    "edgeR_m": "edge_R_m",
    "grade_rad": "grade_rad",
    "isOpen": "isOpen",
    "k_1pm": "k_1pm",
    "posE_m": "posE_m",
    "posN_m": "posN_m",
    "psi_rad": "Psi_rad",
    "s_m": "s_m",
    "UxDes_mps": "Ux_des_mps"
}

def worldfile2pathmsg(fname):
    p = path()
    with open(fname) as f:
        world_dict = yaml.load(f)
        for k in world_dict:
            if isinstance(world_dict[k], str):
                world_dict[k] = map(float, world_dict[k].split(","))
    for k, v in world_dict.items():
        setattr(p, world2path[k], v)
    p.header.frame_id = fname.split(".")[0]
    return p

for fname in glob.glob("*.world"):
    with open(fname.split(".")[0] + ".msg", "w") as f:
        worldfile2pathmsg(fname).serialize(f)
