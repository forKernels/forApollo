// forApollo CUDA bindings — three-state detection
// Copyright The Fantastic Planet - By David Clabaugh

pub extern "c" fn fc_rt_is_stub() c_int;
pub extern "c" fn fc_rt_is_available() c_int;
pub extern "c" fn fc_rt_init() c_int;

pub const CudaAvailability = enum { not_linked, linked_no_gpu, ready };

pub fn detectCuda() CudaAvailability {
    if (fc_rt_is_stub() == 1) return .not_linked;
    if (fc_rt_is_available() == 0) return .linked_no_gpu;
    return .ready;
}

pub fn init() bool { return fc_rt_init() == 0; }
