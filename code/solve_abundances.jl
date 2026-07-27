#!/usr/bin/env julia
# Step 61 — reproduce the catalog from the released-MRT EWs using the VALIDATED
# pipeline routine (Step 30): solve vmic by reduced-EW balance, then [Fe/H]=median Fe I.
# Loads Korg once, loops all 50.  julia 61_solve_from_mrt.jl <manifest> <out_csv>
using Korg, DelimitedFiles, Printf, Statistics, CSV, DataFrames

man = CSV.read(ARGS[1], DataFrame)
ATM_MH = 0.3
out = DataFrame(gid=Int[], inst=String[], src=String[], vmic=Float64[],
                feh_I=Float64[], feh_II=Float64[], ion=Float64[], nI=Int[], scatter=Float64[])

function solve_one(ewcsv, teff, logg)
    data,hdr = readdlm(ewcsv, ',', header=true); cols=vec(hdr)
    wl=Float64.(data[:,findfirst(==("wl_air"),cols)]); sp=string.(data[:,findfirst(==("species"),cols)])
    ew=Float64.(data[:,findfirst(==("ew_mA"),cols)])
    gf=Float64.(data[:,findfirst(==("loggf"),cols)]); xx=Float64.(data[:,findfirst(==("Xexc"),cols)])
    lines=Korg.Line[]; EWs=Float64[]; isFe1=Bool[]
    for i in eachindex(wl)
        species = sp[i]=="Fe I" ? Korg.Species("Fe I") : Korg.Species("Fe II")
        push!(lines, Korg.Line(Korg.air_to_vacuum(wl[i])*1e-8, gf[i], species, xx[i]))
        push!(EWs, ew[i]); push!(isFe1, sp[i]=="Fe I")
    end
    p=sortperm([l.wl for l in lines]); lines=lines[p]; EWs=EWs[p]; isFe1=isFe1[p]
    rew=[log10(EWs[i]*1e-3/Korg.vacuum_to_air(lines[i].wl*1e8)) for i in eachindex(lines)]
    slope(x,y)=cov(x,y)/var(x)
    function eval_vmic(v)
        A_X=Korg.format_A_X(ATM_MH); atm=Korg.interpolate_marcs(teff,logg,A_X)
        A=Korg.Fit.ews_to_abundances(atm, lines, A_X, EWs; vmic=v)
        ok=isfinite.(A); i1=ok.&isFe1; i2=ok.&.!isFe1
        feI=median(A[i1])-7.50; feII=(sum(i2)>0) ? median(A[i2])-7.50 : NaN
        rs=slope(rew[i1],A[i1]); (feI,feII,rs,feI-feII,std(A[i1])/sqrt(sum(i1)),sum(i1))
    end
    lo,hi=0.7,2.8; flo=eval_vmic(lo)[3]; fhi=eval_vmic(hi)[3]
    vbest = sign(flo)==sign(fhi) ? (abs(flo)<abs(fhi) ? lo : hi) : (lo+hi)/2
    if sign(flo)!=sign(fhi)
        for _ in 1:7
            mid=(lo+hi)/2; fm=eval_vmic(mid)[3]
            if sign(fm)==sign(flo); lo=mid; flo=fm; else; hi=mid; fhi=fm; end
            vbest=mid
        end
    end
    feI,feII,rs,id,sc,nI = eval_vmic(vbest)
    (vbest,feI,feII,id,nI,sc)
end

for row in eachrow(man)
    v,fI,fII,id,nI,sc = solve_one(row.ewcsv, row.teff, row.logg)
    push!(out,(row.gid,row.inst,row.src,round(v,digits=2),fI,fII,id,nI,sc))
    @printf("%s %-5s %-6s vmic=%.2f(cat %.2f) feh=%+.3f cat=%+.3f d=%+.3f\n",
            row.gid,row.inst,row.src,v,row.vmic,fI,row.cat_FeH,fI-row.cat_FeH)
end
CSV.write(ARGS[2], out); println("\nwrote ",ARGS[2])
