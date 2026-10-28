module BiochemicalAlgorithmsMetalExt

using BiochemicalAlgorithms
using Metal

# Register Apple Silicon GPUs for the `:gpu` evaluation backend. The actual
# kernels are device-agnostic and live in the base package (ff_gpu.jl); this
# extension only makes the KernelAbstractions device discoverable. Metal does
# not support Float64 in kernels, so Float64 accumulation is disallowed.
function __init__()
    if Metal.functional()
        BiochemicalAlgorithms._register_gpu_device!(Metal.MetalBackend(); supports_f64 = false)
    end
end

end # module
