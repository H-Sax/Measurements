using JLD, CairoMakie
CairoMakie.activate!(type = "svg")
# Load data
begin
    x2 = ones(36)
    x3b = 2*ones(36)
    x3d = 3*ones(36)
    x4b = 4*ones(36)
    x4d = 5*ones(36)
    x4f = 6*ones(36)
    x5b = 7*ones(36)
    x5d = 8*ones(36)
    x5f = 9*ones(36)
    x5h = 10*ones(36)
    x6b = 11*ones(36)
    x7b = 12*ones(36)
    x7d = 13*ones(36)

    y2 =  load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/2.jld")["data"] .+ 1e-16
    y3b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/3b.jld")["data"] .+ 1e-16
    y3d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/3d.jld")["data"] .+ 1e-16
    y4b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/4b.jld")["data"] .+ 1e-16
    y4d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/4d.jld")["data"] .+ 1e-16
    y4f = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/4f.jld")["data"] .+ 1e-16
    y5b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/5b.jld")["data"] .+ 1e-16
    y5d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/5d.jld")["data"] .+ 1e-16
    y5f = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/5f.jld")["data"] .+ 1e-16
    y5h = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/5h.jld")["data"] .+ 1e-16
    y6b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/6b.jld")["data"] .+ 1e-16
    y7b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/7b.jld")["data"] .+ 1e-16
    y7d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/All Levels/7d.jld")["data"] .+ 1e-16
end

f = Figure();
ax = Axis(f[1, 1], xlabel = "Mixed Measurement Sets", ylabel = "Log10 Eigenvalues", yscale = log10, xticks = (1:13, ["2", "3B", "3D", "4B", "4D", "4F", "5B", "5D", "5F", "5H", "6B", "7B", "7D"]))
scatter!(x2,y2, marker = :hline, color = :green)
scatter!(x3b,y3b,marker = :hline, color = :green)
scatter!(x3d,y3d,marker = :hline, color = :green)
scatter!(x4b,y4b,marker = :hline, color = :green)
scatter!(x4d,y4d,marker = :hline, color = :green)
scatter!(x4f,y4f,marker = :hline, color = :green)
scatter!(x5b,y5b,marker = :hline, color = :green)
scatter!(x5d,y5d,marker = :hline, color = :green)
scatter!(x5f,y5f,marker = :hline, color = :green)
scatter!(x5h,y5h,marker = :hline, color = :green)
scatter!(x6b,y6b,marker = :hline, color = :green)
scatter!(x7b,y7b,marker = :hline, color = :green)
scatter!(x7d,y7d,marker = :hline, color = :green)
f