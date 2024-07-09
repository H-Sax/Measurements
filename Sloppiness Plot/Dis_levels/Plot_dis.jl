using JLD, CairoMakie
CairoMakie.activate!(type = "svg")
# Load data
begin
    x1  = 1*ones(36)
    x2a = 2*ones(36)
    x2b = 3*ones(36)
    x2c = 4*ones(36)
    x2d = 5*ones(36)
    x3a = 6*ones(36)
    x3b = 7*ones(36)
    x3c = 8*ones(36)
    x3d = 9*ones(36)
    x3e = 10*ones(36)
    x3f = 11*ones(36)


    y1 =  load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/1.jld")["data"] .+ 1e-16
    y2a = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/2a.jld")["data"] .+ 1e-16
    y2b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/2b.jld")["data"] .+ 1e-16
    y2c = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/2c.jld")["data"] .+ 1e-16
    y2d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/2d.jld")["data"] .+ 1e-16
    y3a = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3a.jld")["data"] .+ 1e-16
    y3b = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3b.jld")["data"] .+ 1e-16
    y3c = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3c.jld")["data"] .+ 1e-16
    y3d = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3d.jld")["data"] .+ 1e-16
    y3e = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3e.jld")["data"] .+ 1e-16
    y3f = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/3f.jld")["data"] .+ 1e-16
end


f = Figure();
ax = Axis(f[1, 1], xlabel = "Discrete Measurement Sets", ylabel = "Log10 Eigenvalues", yscale = log10, xticks = (1:11, ["1", "2A", "2B", "2C", "2D", "3A", "3B", "3C", "3D", "3E", "3F"]))
scatter!(x1 ,y1 , marker = :hline, color = :green)
scatter!(x2a,y2a,marker = :hline, color = :green)
scatter!(x2b,y2b,marker = :hline, color = :green)
scatter!(x2c,y2c,marker = :hline, color = :green)
scatter!(x2d,y2d,marker = :hline, color = :green)
scatter!(x3a,y3a,marker = :hline, color = :green)
scatter!(x3b,y3b,marker = :hline, color = :green)
scatter!(x3c,y3c,marker = :hline, color = :green)
scatter!(x3d,y3d,marker = :hline, color = :green)
scatter!(x3e,y3e,marker = :hline, color = :green)
scatter!(x3f,y3f,marker = :hline, color = :green)
f