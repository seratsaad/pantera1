#!/usr/bin/env julia
# Step 101 -- Build the public release table of iron abundances at four
# equivalent-width caps (full, 140, 120, 100 mA) for all 62 spectra.
#
# Self-contained: reads only the released tables, so the result is reproducible
# from the published data alone.
#   pantera_letter/table_ews_mrt.csv      the released equivalent widths
#   pantera_letter/table_results_mrt.csv  adopted Teff, logg, published FeH
#
# For each spectrum and each cap we redo the full abundance step: vmic from the
# reduced-equivalent-width balance over the retained Fe I lines, then
# [Fe/H] = median Fe I abundance - 7.50. Same machinery as Step 30 and Step 80.
#
# Out: results/release/pantera_feh_ewcaps.csv   (per spectrum)
#      results/release/pantera_feh_ewcaps_star.csv (per star, mean over spectra)
using Korg, CSV, DataFrames, Statistics, Printf

ROOT = "/Users/saad.104/Downloads/pantera"
OUTD = "$ROOT/results/release"
mkpath(OUTD)

ews = CSV.read("$ROOT/pantera_letter/table_ews_mrt.csv", DataFrame)
res = CSV.read("$ROOT/pantera_letter/table_results_mrt.csv", DataFrame)
CAPS = [("full", Inf), ("w140", 140.0), ("w120", 120.0), ("w100", 100.0)]

out = DataFrame(gaia_dr3=Int[], instrument=String[], cap=String[],
                teff=Float64[], logg=Float64[], vmic=Float64[],
                FeH=Float64[], FeII=Float64[], nFeI=Int[], nFeII=Int[],
                railed=Int[])

outfile = "$OUTD/pantera_feh_ewcaps.csv"
println("solving $(nrow(res)) spectra x $(length(CAPS)) caps ...")

for r in eachrow(res)
    sub = ews[(ews.gaia_dr3 .== r.gaia_dr3) .& (ews.instrument .== r.instrument), :]
    nrow(sub) == 0 && (println("no EWs for $(r.gaia_dr3) $(r.instrument)"); continue)
    mh0 = clamp(r.FeH_dex, 0.0, 0.45)
    A_X = Korg.format_A_X(mh0)
    atm = Korg.interpolate_marcs(r.teff_K, r.logg_dex, A_X)   # same for every cap

    for (tag, cap) in CAPS
        s = sub[sub.ew_mA .<= cap, :]
        nfe1 = sum(s.species .== "Fe I")
        nfe1 < 15 && (println("  skip $(r.gaia_dr3) $tag, only $nfe1 Fe I lines"); continue)

        lines = Korg.Line[]; EWs = Float64[]; isFe1 = Bool[]
        for l in eachrow(s)
            sp = l.species == "Fe I" ? Korg.Species("Fe I") : Korg.Species("Fe II")
            push!(lines, Korg.Line(Korg.air_to_vacuum(l.wl_air_AA)*1e-8, l.loggf, sp, l.chi_eV))
            push!(EWs, l.ew_mA); push!(isFe1, l.species == "Fe I")
        end
        p = sortperm([l.wl for l in lines])
        lines, EWs, isFe1 = lines[p], EWs[p], isFe1[p]
        rew = [log10(EWs[i]*1e-3/Korg.vacuum_to_air(lines[i].wl*1e8)) for i in eachindex(lines)]
        slope(x, y) = cov(x, y)/var(x)

        function ev(v)
            A = Korg.Fit.ews_to_abundances(atm, lines, A_X, EWs; vmic=v)
            ok = isfinite.(A); i1 = ok .& isFe1; i2 = ok .& .!isFe1
            fe1 = median(A[i1]) - 7.50
            fe2 = sum(i2) > 0 ? median(A[i2]) - 7.50 : NaN
            fe1, slope(rew[i1], A[i1]), sum(i1), fe2, sum(i2)
        end

        lo, hi = 0.7, 2.8
        flo = ev(lo)[2]; fhi = ev(hi)[2]
        railed = 0; v = (lo+hi)/2
        if sign(flo) == sign(fhi)
            v = abs(flo) < abs(fhi) ? lo : hi; railed = 1
        else
            for _ in 1:7
                fm = ev(v)[2]
                if sign(fm) == sign(flo); lo = v; flo = fm; else; hi = v; end
                v = (lo+hi)/2
            end
        end
        fe1, _, n1, fe2, n2 = ev(v)
        push!(out, (r.gaia_dr3, r.instrument, tag, r.teff_K, r.logg_dex,
                    round(v, digits=3), round(fe1, digits=4),
                    isnan(fe2) ? NaN : round(fe2, digits=4), n1, n2, railed))
        @printf("%d %-5s %-4s vmic=%.2f feh=%+.3f nFeI=%d%s\n",
                r.gaia_dr3, r.instrument, tag, v, fe1, n1, railed == 1 ? "  RAILED" : "")
        CSV.write(outfile, out)      # incremental, so partial results survive
    end
end

CSV.write(outfile, out)

# per-star means over instruments
star = combine(groupby(out, [:gaia_dr3, :cap]),
               :FeH => mean => :FeH, :vmic => mean => :vmic,
               :nFeI => mean => :nFeI, :railed => maximum => :railed_any,
               nrow => :n_spectra)
star.FeH = round.(star.FeH, digits=4); star.vmic = round.(star.vmic, digits=3)
CSV.write("$OUTD/pantera_feh_ewcaps_star.csv", star)

println("\nrows: $(nrow(out))  stars: $(length(unique(out.gaia_dr3)))")
for (tag, _) in CAPS
    s = out[out.cap .== tag, :]
    println("  $tag: $(nrow(s)) spectra, median FeH $(round(median(s.FeH), digits=3)), railed $(sum(s.railed))")
end
println("wrote $outfile and the per-star version")
