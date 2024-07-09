using JLD, CairoMakie
CairoMakie.activate!(type = "svg")
# Load data
begin
    x1a = 1*ones(36)
    x1b = 2*ones(36)
    x1c = 3*ones(36)
    x1d = 4*ones(36)
    x1e = 5*ones(36)
    x1f = 6*ones(36)
    x2a = 7*ones(36)
    x2b = 8*ones(36)
    x2c = 9*ones(36)
    x2d = 10*ones(36)
    x3a = 11*ones(36)
    x3b = 12*ones(36)
    x4a = 13*ones(36)
    x4b = 14*ones(36)
    x4c = 15*ones(36)
    x4d = 16*ones(36)

    y1a=  load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1a.jld")["data"] .+ 1e-16
    y1b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1b.jld")["data"] .+ 1e-16
    y1c = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1c.jld")["data"] .+ 1e-16
    y1d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1d.jld")["data"] .+ 1e-16
    y1e = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1e.jld")["data"] .+ 1e-16
    y1f = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/1f.jld")["data"] .+ 1e-16
    y2a = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/2a.jld")["data"] .+ 1e-16
    y2b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/2b.jld")["data"] .+ 1e-16
    y2c = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/2c.jld")["data"] .+ 1e-16
    y2d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/2d.jld")["data"] .+ 1e-16
    y3a = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/3a.jld")["data"] .+ 1e-16
    y3b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/3b.jld")["data"] .+ 1e-16
    y4a = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/4a.jld")["data"] .+ 1e-16
    y4b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/4b.jld")["data"] .+ 1e-16
    y4c = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/4c.jld")["data"] .+ 1e-16
    y4d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Con_Levels/4d.jld")["data"] .+ 1e-16
end


f = Figure();
ax = Axis(f[1, 1], xlabel = "Continous Measurement Sets", ylabel = "Log10 Eigenvalues", yscale = log10, xticks = (1:16, ["1A", "1B", "1C", "1D", "1E", "1F", "2A", "2B", "2C", "2D", "3A", "3B", "4A","4B","4C","4D"]))
scatter!(x1a,y1a, marker = :hline, color = :green)
scatter!(x1b,y1b,marker = :hline, color = :green)
scatter!(x1c,y1c,marker = :hline, color = :green)
scatter!(x1d,y1d,marker = :hline, color = :green)
scatter!(x1e,y1e,marker = :hline, color = :green)
scatter!(x1f,y1f,marker = :hline, color = :green)
scatter!(x2a,y2a,marker = :hline, color = :green)
scatter!(x2b,y2b,marker = :hline, color = :green)
scatter!(x2c,y2c,marker = :hline, color = :green)
scatter!(x2d,y2d,marker = :hline, color = :green)
scatter!(x3a,y3a,marker = :hline, color = :green)
scatter!(x3b,y3b,marker = :hline, color = :green)
scatter!(x4a,y4a,marker = :hline, color = :green)
scatter!(x4b,y4b,marker = :hline, color = :green)
scatter!(x4c,y4c,marker = :hline, color = :green)
scatter!(x4d,y4d,marker = :hline, color = :green)
f