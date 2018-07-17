using Pigeon
using PyCall

function load_path_msg(fname, msg=PyCall.PyObject(Pigeon.osprey.msg.path()))
    open(fname) do f
        Pigeon.osprey.msg.path[:deserialize](msg, read(f))
        msg[:header][:frame_id] = String(msg[:header][:frame_id])
        Pigeon.osprey.msg.path(msg)
    end
end

function load_all_path_msgs(dir=joinpath(Pkg.dir("Pigeon"), "test", "path"))
    Dict(fname => load_path_msg(joinpath(dir, fname)) for fname in readdir(dir) if endswith(fname, ".msg"))
end

@eval(Pigeon, const TEST_PATHS = $(load_all_path_msgs()))
