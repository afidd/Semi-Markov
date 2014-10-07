using DataFrames
using Gadfly

if length(ARGS)<3
    println("julia comparison.jl individual_cnt R0 run_cnt")
    exit()
end
cnt=int(ARGS[1])
r0=float(ARGS[2])
runs=int(ARGS[3])

function analytical_distribution(cnt, r0)
    expected=zeros(Float64, cnt)
    open(`python Ball_Nasell.py -N $(cnt-1) -a 1 -r $(r0)`, "r") do io
        two_numbers=r"^([0-9]+), ([0-9.]+)"
        for l in eachline(io)
            m=match(two_numbers, chomp(l))
            if m!=nothing
                expected[int(m.captures[1])+1]=float(m.captures[2])
            end
        end
    end
    expected
end

function simulated_distribution(cnt, r0, runs)
    histogram=zeros(Int,cnt)
    b=r0/cnt # Translating from well-mixed to fully-expanded beta.
    open(`../sirmixed --run $(runs) -j 8 --loglevel=error -s $(cnt) --beta=$(b) --gamma=1`, "r") do f
        line_number=0
        for l in eachline(f)
            l=chop(l)
            if length(l)<1
                println(line_number, ":",l,":")
            else
                try
                    value=int(l)
                    histogram[value]=histogram[value]+1
                catch
                    println("could not interpret :", l, ":")
                end
            end
            line_number+=1
        end
    end
    histogram
end

function make_pdf(analytical, simulated, pdfname)
    x=reshape(1:cnt, cnt)
    assert(length(x)==length(simulated))
    df=DataFrame(Recovered=x, Count=simulated/sum(simulated), Expected=analytical)
    myplot=plot(df, layer(x="Recovered", y="Count", Geom.point,
            Theme(default_color=color("blue"))),
        layer(x="Recovered", y="Expected", Geom.point,
            Theme(default_color=color("orange"))),
        Guide.xlabel("Infected"), Guide.ylabel("PDF"),
        Guide.title("SIR Found versus Analytic Result"))
    draw(PDF("$(pdfname).pdf", 20cm, 15cm), myplot)
end

simulated=simulated_distribution(cnt, r0, runs)
analytical=analytical_distribution(cnt, r0)
make_pdf(analytical, simulated, "CheckOnSIR")

