
"""
plot!(p_main, 
    [0, 0.01], [1.6, 1.6], 
    subplot = 1, 
    line = (:dash, :black, 1.5), 
    label = ""
)

plot!(p_main, 
    [0.01, 0.01], [4/3, 1.6], 
    subplot = 1, 
    line = (:dash, :black, 1.5), 
    label = ""
)

using Plots: bbox  # to access bbox()

# 1. Define the inset within the main plot call
plot!(p_main, 
    Q_y[1:11], Q_x[1:11],
    inset = bbox(0.54, 0.38, 0.38, 0.40), # Position and size
    subplot = 2,                          # Directs the following commands to the inset
    xlims = (0.0, 0.01),
    ylims = (1.333, 1.601),
    xticks = 0:0.002:0.01,                # Adjusted for the smaller scale
    #yticks = 0:0.002:0.01,                # Adjusted for the smaller scale
    framestyle = :box,
    legend = :none,
    xtickfont = font(12),                 # Scaled down for readability in small box
    ytickfont = font(12),
    guidefont = font(14),
    xlabel = "", 
    ylabel = ""
)
# Connection from Top-Right of region (0.01, 1.6) to Bottom-Left of Inset
plot!(p_main, [0.01, 0.059], [1.6, 1.717], line=(:dot, :gray, 1), subplot=1, label="")

plot!(p_main, [0, 0.0315], [1.6, 1.717], line=(:dot, :gray, 1), subplot=1, label="")

# Connection from Bottom-Right of region (0.01, 4/3) to Bottom-Left of Inset
plot!(p_main, [0.01, 0.059], [4/3, 1.386], line=(:dot, :gray, 1), subplot=1, label="")

plot!(p_main, [0, 0.0315], [4/3, 1.386], line=(:dot, :gray, 1), subplot=1, label="")

display(p_main)
savefig("Ent_imprec_pauli.pdf")
"""

using LaTeXStrings
using Plots
using DelimitedFiles  
using FileIO, JLD2    

Q_data = readdlm("ent_imprec_data.txt")      # change delimiter if needed
Q_y = Q_data[:, 1]
Q_x = Q_data[:, 2]

gr(size = (1000, 500), linewidth = 4)
# Plot the classical data:
p_main = plot(
    Q_y,  
    Q_x,  
    legend = :none, 
    framestyle = :box,
    xlims    = (0, 0.063),
    ylims    = (4/3, 2.01),
    xticks = 0:0.01:0.064,
    left_margin=5.5Plots.mm,
    right_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
    xlabel = L"ε", 
    ylabel = L"W",
    color=colorant"#0B299E",
    xtickfont=font(20), 
    ytickfont=font(20), 
    guidefont=font(26), 
    legendfont=font(22)
)


hline!([2], linestyle=:dash, color=:gray, label=false)

savefig("Ent_imprec_pauli.pdf")