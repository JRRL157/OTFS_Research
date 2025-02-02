//OTFS Frame parameters

N = 16;
M = 64;

Fn = dftmtx(N);
Fn=Fn/norm(Fn);

delta_f = 15e3;

T=1/delta_f;

fc = 4e9;

c = 299792458;

delay_resolution = 1/(M*delta_f);
doppler_resolution = 1/(N*T);


// Generating OTFS frame
mod_size = 4;

N_syms_per_frame = N*M;

N_bits_per_frame = N*M*log2(mod_size);

