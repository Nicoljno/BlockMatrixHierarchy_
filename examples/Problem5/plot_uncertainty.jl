using LaTeXStrings
using Plots
using DelimitedFiles  
using FileIO, JLD2    

Q_data = readdlm("uncertainty_N5_d2.txt")      # change delimiter if needed
Q_y = Q_data[:, 1]
Q_x = Q_data[:, 2]

Q3_data = readdlm("uncertainty_N3_d2.txt")      # change delimiter if needed
Q3_y = Q3_data[:, 1]
Q3_x = Q3_data[:, 2]

Q4_data = readdlm("uncertainty_N4_d2.txt")      # change delimiter if needed
Q4_y = Q4_data[:, 1]
Q4_x = Q4_data[:, 2]

gr(size = (1000, 500), linewidth = 4)
# Plot the classical data:
plot(
    Q_y,  
    Q_x,
    framestyle = :box,
    xlims    = (0, 1),
    xticks = 0:0.1:1,
    left_margin=5.5Plots.mm,
    right_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
    label = "N=5",
    xlabel = L"ε", 
    ylabel = L"β",
    xtickfont=font(20), 
    ytickfont=font(20), 
    guidefont=font(26), 
    legendfont=font(20)
)


plot!(
    Q4_y,  
    Q4_x,
    label = "N=4"
)
plot!(
    Q3_y,  
    Q3_x,
    label = "N=3"
)
savefig("plot_uncertainty.pdf")