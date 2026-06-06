module BiochemicalAlgorithmsCUDAExt

using BiochemicalAlgorithms
using CUDA

# Register NVIDIA GPUs for the `:gpu` evaluation backend. The kernels are
# device-agnostic and live in the base package (ff_gpu.jl); this extension only
# makes the KernelAbstractions device discoverable. CUDA supports Float64.
function __init__()
    if CUDA.functional()
        BiochemicalAlgorithms._register_gpu_device!(CUDA.CUDABackend(); supports_f64 = true)
    end
end

end # module
