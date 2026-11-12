#!/usr/bin/env python3
"""RAINBOR: sparse 12-billion-address interference and exact freeze/play jump test.

The Fraunhofer screen calculation follows RAINBOR-DOUBLE-SLIT.html:
 single-slit sinc^2 * normalized N-slit interference for 1,2,3 slits.
 A 409-channel synthetic spectrum is *not* 409 archived GGUF specimens.
 The 27-site linear dynamical state is imported from the source-pinned
 matrix-law essentials code (the NTT solver and exact modular field).
 Source density, temperature cutoff, wind, radiation and state-to-phase
 coupling are explicitly chosen exploratory optical rules, NOT gravity,
 Hawking radiation, quantum particles, or a Standard Model calculation.
 Only 243 of 12,000,000,000 possible virtual screen indices are sampled.
 No 12B-sized vector/matrix is ever constructed. No coupling-r9 code.
"""
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
from hashlib import sha256
from math import pi, sin, cos, sqrt, log1p, exp
import argparse
import tracemalloc

VIRTUAL_PIXELS=12_000_000_000
SAMPLES=243
P=1_000_081
LAM0=18.0
SEP=90.0
WIDTH=8.0
K40=3**40
ASOLARIA_SURROGATE_CHANNELS=409

def load_jump(path):
    spec=spec_from_file_location("rainbor_matrix_law_source_pinned",str(path))
    assert spec and spec.loader
    mod=module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.P==P
    return mod

def pixels():
    out=[i*(VIRTUAL_PIXELS-1)//(SAMPLES-1) for i in range(SAMPLES)]
    assert out[0]==0 and out[-1]==VIRTUAL_PIXELS-1 and len(set(out))==SAMPLES
    return out

def u_from_pixel(idx):
    assert 0<=idx<VIRTUAL_PIXELS
    return 2.0*idx/(VIRTUAL_PIXELS-1)-1.0

def sinc(x):
    return 1.0 if abs(x)<1e-14 else sin(x)/x

def original_screen_i(u,slits,lam=LAM0):
    """Direct translation of RAINBOR-DOUBLE-SLIT.html function screenI()."""
    assert slits in (1,2,3)
    sintheta=u*0.5
    beta=pi*WIDTH*sintheta/lam
    gamma=pi*SEP*sintheta/lam
    single=sinc(beta)**2
    multi=1.0
    if slits>1:
        sg=sin(gamma)
        multi=1.0 if abs(sg)<1e-6 else (sin(slits*gamma)/sg)**2/(slits*slits)
    return single*multi

def state_phase(state,slit):
    # Slit k taps a different one of 27 exact modular signal channels.
    # This is an experimental coupling, not derived from Asolaria photons.
    return 0.24*(state[(9*slit+3)%27]-state[(9*slit+8)%27])/P

def effective_source(rho,temp_over_cap,gate):
    assert rho>=0 and temp_over_cap>=0
    cutoff=max(0.0,1.0-temp_over_cap) if gate else 1.0
    return 0.06*rho*cutoff

def intensity(u,slits,lam,state,rho=0.0,temp_over_cap=0.0,gate=True,
              wind=0.0,radiation=False,drive=True):
    assert slits in (1,2,3) and lam>0
    beta=pi*WIDTH*(u*0.5)/lam
    gamma=pi*SEP*(u*0.5)/lam
    coords=([0.0] if slits==1 else [-0.5,0.5] if slits==2 else [-1.0,0.0,1.0])
    strength=effective_source(rho,temp_over_cap,gate)
    re=im=0.0
    for s,j in enumerate(coords):
        # Softened, off-axis logarithmic lens: dimensionless phase *proxy*.
        lens=strength*(LAM0/lam)*log1p(
            1.0 / sqrt(0.06**2+(u-0.17)**2+(j-0.19)**2)
        ) if strength else 0.0
        modulation=state_phase(state,s) if drive else 0.0
        phase=2.0*j*gamma + lens + 0.16*wind*j*u + modulation
        re+=cos(phase)
        im+=sin(phase)
    val=sinc(beta)**2*(re*re+im*im)/(slits*slits)
    if radiation:
        # Additive *uncorrelated* synthetic source halo, not Hawking emission.
        val+=0.025*temp_over_cap*rho/(rho+12.0)*exp(-((u-0.17)/0.28)**2)
    return val

def state_hash(jump,state):
    return jump.state_sha(state)

def run_jumps(jump,out):
    ring=jump.Ring(27,((0,1),(1,1),(2,-1)))
    ops=jump.Ops()
    lam=jump.eigenvalues(ring,ops)
    assert len(lam)==27 and all(x!=0 for x in lam)
    x0=jump.start_state(27)
    assert len(x0)==27
    known={3:"7d05096afc24b08c",6:"cc3581ab163de3a7",
           9:"ee4b40fb753713ab",10:"81323c66a0e8a92f",
           20:"8ef1ccc86cd3ff0c",40:"177066f133a5c9ef"}
    ladder=jump.cube_ladder(jump.circulant(ring),40,jump.Ops())
    saved={"t0":x0}
    for k in (3,6,9,10,20,40):
        n=3**k
        state=jump.ntt_jump(ring,lam,x0,n,jump.Ops())
        h=state_hash(jump,state)
        assert h.startswith(known[k]),("independent matrix-law published hash",k,h)
        cubed=jump.matvec(ladder[k][0],x0,jump.Ops())
        assert cubed==state,("NTT cube mismatch",k)
        if k<=6:
            direct=x0
            for _ in range(n):
                direct=jump.step(ring,direct,jump.Ops())
            assert direct==state,("direct mismatch",k)
            direct_verdict="PASS"
        else:direct_verdict="NOT_RUN"
        if k in (6,10,40):saved["3pow"+str(k)]=state
        out.append("JUMP|method=NTT_AND_CUBE|k=%d|steps=3^%d|state_sha256=%s|cube_eq_ntt=PASS|direct=%s|integer=EXACT|json=0"%
                   (k,k,h,direct_verdict))
    t40=K40
    x40=saved["3pow40"]
    fwd=jump.step(ring,x40,jump.Ops())
    assert fwd==jump.ntt_jump(ring,lam,x0,t40+1,jump.Ops())
    fwd3=x40
    for _ in range(3):fwd3=jump.step(ring,fwd3,jump.Ops())
    assert fwd3==jump.ntt_jump(ring,lam,x0,t40+3,jump.Ops())
    prev=jump.ntt_unstep(ring,lam,x40,jump.Ops())
    assert prev is not None
    assert jump.step(ring,prev,jump.Ops())==x40
    assert prev==jump.ntt_jump(ring,lam,x0,t40-1,jump.Ops())
    saved["3pow40+1"]=fwd
    saved["3pow40+3"]=fwd3
    out.append("FREEZE_PLAY|at=3^40|forward1=PASS|forward3=PASS|backward1=PASS|state_sha256=%s|integer=EXACT|json=0"%
               state_hash(jump,x40))
    t_binary=2**40
    binary=jump.ntt_jump(ring,lam,x0,t_binary,jump.Ops())
    saved["2pow40"]=binary
    out.append("BINARY_JUMP|steps=2^40|state_sha256=%s|integer=EXACT|json=0"%
               state_hash(jump,binary))
    off=2**20+3
    from_freeze=jump.ntt_jump(ring,lam,x40,off,jump.Ops())
    straight=jump.ntt_jump(ring,lam,x0,t40+off,jump.Ops())
    assert from_freeze==straight
    saved["mixed"]=straight
    out.append("MIXED_JUMP|time=3^40+2^20+3|frozen_then_binary_eq_direct_NTT=PASS|state_sha256=%s|integer=EXACT|json=0"%
               state_hash(jump,straight))
    # NEGATIVE CONTROL: a local density kick applied at every step does NOT
    # commute with the original pure M power. Jumping M alone is invalid.
    ordinary=jump.ntt_jump(ring,lam,x0,27,jump.Ops())
    disturbed=x0[:]
    for _ in range(27):
        disturbed=jump.step(ring,disturbed,jump.Ops())
        disturbed[0]=(disturbed[0]+7)%P
    assert disturbed!=ordinary
    out.append("NONLINEAR_OR_FORCED_CONTROL|repeated_local_source_KICK=7|steps=27|pure_linear_jump_eq_forced=NO|correctly_refused=PASS|json=0")
    # Negative cache/receipt check: any changed state needs a different digest.
    tampered=x40[:];tampered[0]=(tampered[0]+1)%P
    assert state_hash(jump,tampered)!=state_hash(jump,x40)
    out.append("HASH_MUTATION|change_one_component=1|sha256_changes=PASS|hash_is_not_force=1|json=0")
    return saved

def run_baseline(out):
    error=0.0
    for slits in (1,2,3):
        for z in range(2001):
            u=-1.0+2.0*z/2000
            a=original_screen_i(u,slits)
            b=intensity(u,slits,LAM0,[0]*27,drive=False)
            error=max(error,abs(a-b))
    assert error<1e-12,("RAINBOR HTML baseline mismatch",error)
    out.append("LEGACY_FRAUNHOFER|slits=1,2,3|angles=2001_each|maximum_abs_error=%.12g|gate=PASS|json=0"%error)

def summarize(values):
    ave=sum(values)/len(values)
    npeak=sum(values[i]>values[i-1] and values[i]>values[i+1] for i in range(1,len(values)-1))
    center=values[len(values)//2]
    skew=sum(values[(SAMPLES//2)+i]-values[(SAMPLES//2)-i] for i in range(1,SAMPLES//2))
    return ave,center,npeak,skew/len(values)

def rms(a,b):
    assert len(a)==len(b) and len(a)>0
    return sqrt(sum((x-y)**2 for x,y in zip(a,b))/len(a))

def render(state,slits,spectrum,case,wind=0.0,radiation=False):
    rho,temp,gate=case
    addrs=pixels()
    lambdas=[LAM0] if spectrum=="mono" else [8.0+32.0*i/408 for i in range(409)]
    out=[]
    for addr in addrs:
        u=u_from_pixel(addr)
        total=0.0
        for lam in lambdas:
            total+=intensity(u,slits,lam,state,rho,temp,gate,wind,radiation,True)
        out.append(total/len(lambdas))
    return out

def experiment(saved,out):
    cases={
        "vacuum":(0.0,0.0,True),
        "ordinary":(1.0,0.10,True),
        "neutron_like":(12.0,0.25,True),
        "blackhole_like":(80.0,0.65,True),
        "supercritical_cutoff":(250.0,1.20,True),
        "supercritical_no_cutoff":(250.0,1.20,False),
    }
    cases_order=tuple(cases)
    assert effective_source(*cases["supercritical_cutoff"])==0
    selected=[("initial","t0"),("ternary","3pow40"),("binary","2pow40"),("mixed","mixed")]
    sample_rows=0
    strongest=0.0
    min_lens_difference=None
    for stage,key in selected:
        state=saved[key]
        for slits in (1,2,3):
            for spectrum in ("mono","rainbow409"):
                baseline=render(state,slits,spectrum,cases["vacuum"])
                for name in cases_order:
                    vals=render(state,slits,spectrum,cases[name])
                    error=rms(vals,baseline)
                    av,centre,peaks,skew=summarize(vals)
                    if name=="supercritical_cutoff":
                        assert error<1e-13,("temperature gate did not erase synthetic lens",error)
                    if name=="neutron_like" and slits>1:
                        strongest=max(strongest,error)
                        min_lens_difference=error if min_lens_difference is None else min(min_lens_difference,error)
                    if name=="blackhole_like" and slits==1:
                        assert error<1e-12,("one slit phase-only null failed",error)
                    out.append("RAINBON|time=%s|slits=%d|spectrum=%s|source=%s|sampled=%d|mean=%.9f|centre=%.9f|rms_vs_vacuum=%.9f|sampled_peaks=%d|asymmetry=%.9f|cutoff=%s|json=0"%
                               (stage,slits,spectrum,name,SAMPLES,av,centre,error,peaks,skew,"ON" if cases[name][2] else "OFF"))
                    sample_rows+=1
                # Wind and noncoherent additive background are separate
                # controls evaluated at actual neutron-like source strength.
                ncase=cases["neutron_like"]
                wplus=render(state,slits,spectrum,ncase,wind=+1.0)
                wminus=render(state,slits,spectrum,ncase,wind=-1.0)
                halo=render(state,slits,spectrum,ncase,radiation=True)
                nohalo=render(state,slits,spectrum,ncase,radiation=False)
                winddiff=rms(wplus,wminus)
                radiff=rms(halo,nohalo)
                if slits>1:
                    assert winddiff>1e-8,("wind sign didn't change footprint",stage,slits,spectrum)
                else:
                    assert winddiff<1e-12,("single-slit wind phase null",winddiff)
                assert radiff>0
                out.append("EXTERNAL_CONTROL|time=%s|slits=%d|spectrum=%s|wind_reversal_rms=%.9f|radiation_halo_rms=%.9f|wind_is_assumed=1|radiation_is_assumed=1|json=0"%
                           (stage,slits,spectrum,winddiff,radiff))
    assert sample_rows==4*3*2*len(cases)
    assert strongest>1e-6 and min_lens_difference>1e-8,("density source had no interference response",strongest,min_lens_difference)
    out.append("EXPERIMENT_PASS|scenarios=%d|screens=%d|meta_pixel_addresses=%d|sampled_addresses_per_screen=%d|matrix_jump_integer=EXACT|density_response_nonzero=PASS|temperature_cutoff_null=PASS|one_slit_null=PASS|wind_and_radiation_controls=PASS|json=0"%
               (len(cases),sample_rows,VIRTUAL_PIXELS,SAMPLES))
    return sample_rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--matrix-source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    tracemalloc.start()
    result=[
        "META|schema=RAINBON-METAPIXEL-12B-JUMP-DENSITY-V1|simulation=OPTICAL_MODEL|virtual_addresses=%d|sampled=%d|colors_rainbow=409|r9=HOLD|json=0"%
        (VIRTUAL_PIXELS,SAMPLES)
    ]
    # Source pin goes into the reproducible receipt, before calculations.
    b=args.matrix_source.read_bytes()
    result.append("SOURCE|matrix_law_jump_file_sha256=%s|source_is_unmodified=YES|json=0"%sha256(b).hexdigest())
    jump=load_jump(args.matrix_source)
    saved=run_jumps(jump,result)
    run_baseline(result)
    screens=experiment(saved,result)
    _,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak<30_000_000,("not sparse enough",peak)
    result.append("RESOURCE|virtual_addresses=%d|screens=%d|peak_python_alloc_bytes=%d|dense_12b_allocation=NO|json=0"%
                  (VIRTUAL_PIXELS,screens,peak))
    result.append("LIMIT|rainbons_are_not_actual_photons|inferred_quarks_gluons_glueballs_protons=NO|gravity_temperature_lens_RADIATION_are_assumed_rules|hash_not_gravity|12b_independent_particles=NO|r9_HOLD=YES|json=0")
    result.append("RAINBON12B_VERIFY|status=PASS|matrix_jumps=EXACT|legacy_optics=PASS|density_temperature=ASSUMED|virtual_screen=SPARSE|json=0")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text("\n".join(result)+"\n",encoding="ascii")
    args.output.with_suffix(args.output.suffix+".sha256").write_text(
        sha256(args.output.read_bytes()).hexdigest()+"  "+args.output.name+"\n",encoding="ascii")
    print("\n".join(r for r in result if r.startswith(("JUMP|","BINARY_JUMP|","MIXED_JUMP|","FREEZE_PLAY|","LEGACY_FRAUNHOFER|","EXPERIMENT_PASS|","RESOURCE|","RAINBON12B_VERIFY|"))))

if __name__=="__main__":
    main()
