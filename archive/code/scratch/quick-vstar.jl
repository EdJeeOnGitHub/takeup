using Plots
using StatsBase
using Distributions
using QuadGK
using ForwardDiff
using NLsolve

# Define parameters
μ1, σ1 = -1.0, 1.0  # Mean and standard deviation for the first Gaussian component
μ2, σ2 = 1.0, 2.0  # Mean and standard deviation for the second Gaussian component
π = 0.5  # Mixing proportion

# Function to sample from the mixture model
function sample_gmm(n)
    samples = zeros(n)
    for i in 1:n
        if rand() < π
            samples[i] = rand(Normal(μ1, σ1))  # Sample from the first Gaussian
        else
            samples[i] = rand(Normal(μ2, σ2))  # Sample from the second Gaussian
        end
    end
    return samples
end

# Example: Sample 1000 points from the GMM
n = 10_000
samples = sample_gmm(n)



# Parameters of the GMM
μ1, σ1, π1 = 0.4, 0.2, 0.75
μ2, σ2, π2 = -2.0, 0.2, 0.25

model_params = (
    μ1 = μ1, 
    σ1 = σ1, 
    π1 = π1, 
    μ2 = μ2, 
    σ2 = σ2, 
    π2 = π2
    )
# Define the PDF of the GMM
function gmm_pdf(v, params)
    (μ1, σ1, π1, μ2, σ2, π2) = params
    return π1 * pdf(Normal(μ1, σ1), v) + π2 * pdf(Normal(μ2, σ2), v)
end


function gmm_cdf(v, params)
    (μ1, σ1, π1, μ2, σ2, π2) = params
    return π1 * cdf(Normal(μ1, σ1), v) + π2 * cdf(Normal(μ2, σ2), v)
end


function normal_pdf(v, params)
    (μ1, σ1, π1, μ2, σ2, π2) = params
    return pdf(Normal(μ1, σ1), v) 
end

function normal_cdf(v, params)
    (μ1, σ1, π1, μ2, σ2, π2) = params
    return cdf(Normal(μ1, σ1), v) 
end





pdf(Normal(0.2, 1.0), 1.0)

# Define the conditional expectations functions
function MPlus(x, params, pdf)
    integral, _ = quadgk(v -> v * pdf(v, params), x, Inf)
    probability, _ = quadgk(y -> pdf(y, params), x, Inf)
    return integral / probability
end


function MMinus(x, params, pdf)
    integral, _ = quadgk(v -> v * pdf(v, params), -Inf, x)
    probability, _ = quadgk(y -> pdf(y, params), -Inf, x)
    return integral / probability
end

# Example usage
x = 2.0
println("E[v | v > $x]: ", MPlus(x, model_params, gmm_pdf))
println("E[v | v < $x]: ", MMinus(x, model_params, gmm_pdf))


function Delta(x, params, pdf)
    return MPlus(x, params, pdf) - MMinus(x, params, pdf)
end

# 1 structural_clus… -0.412 
# 2 obs_cluster_mu_…  0.0876
# 3 total_error_sd    1.23  
# 4 u_sd              0.623 

x = -3:0.01:3
μ = 0.0876
b = -0.412
c = 0
plot(x, gmm_pdf.(x, Ref(model_params)), label="Mixture Model PDF", xlabel="\$V^*\$", ylabel="Density", title="\$V^*\$ PDFs")
plot!(x, normal_pdf.(x, Ref(model_params)), label="Gaussian PDF")

savefig(
    "temp-data/gauss-mix-pdf.pdf"
)

plot(x, Delta.(x, Ref(model_params), gmm_pdf), label="Mixture \$\\Delta\$", xlabel="\$v^*\$", ylabel="Reputational Return", title="\$\\Delta(v^*)\$")
plot!(x, Delta.(x, Ref(model_params), normal_pdf), label="Gaussian \$\\Delta\$")

savefig(
    "temp-data/Delta-vstar-mix-pdf.pdf"
)

# Define DeltaPrime(x), the derivative of Delta with respect to x
function DeltaPrime(x, params, pdf; h=1e-5)
    return (Delta(x + h, params, pdf) - Delta(x - h, params, pdf)) / (2h)
end

function FindVStar(params, fp_vector, pdf)
    b, c, μ = fp_vector
    f(x) = x[1]  - c + b + μ * Delta(x[1], params, pdf)
    sol = nlsolve(f, [0.0])
    return sol
end

fp_vec = [b, c, μ]

sol = FindVStar(model_params, fp_vec, gmm_pdf)
bs = -2.5:0.01:2.5
gmm_vstars = map(bs) do b
    FindVStar(model_params, [b, c, μ], gmm_pdf).zero[1]
end
normal_vstars = map(bs) do b
    FindVStar(model_params, [b, c, μ], normal_pdf).zero[1]
end



plot(
    bs, 
    gmm_vstars, 
    label="Bimodal \$V^*(b)\$", 
    xlabel="\$b\$", 
    ylabel="\$V^*(b)\$", 
    title="\$V^*\$ as \$b\$ Changes")
plot!(
    bs, 
    normal_vstars, 
    label="Gaussian \$V^*(b)\$")
savefig(
    "temp-data/Vstar-b.pdf"
)


plot(
    bs,
    1 .- gmm_cdf.(gmm_vstars, Ref(model_params)),
    label = "Bimodal",
    xlabel = "\$b\$",
    ylabel = "Take-up",
    title = "Take-up as \$b\$ Changes"
)
plot!(
    bs,
    1 .- normal_cdf.(normal_vstars, Ref(model_params)),
    label = "Normal",
    xlabel = "\$b\$",
    ylabel = "Take-up"
)

savefig("temp-data/b-changes.pdf")
μs = 0:0.01:1
μ_gmm_vstars = map(μs) do μ
    FindVStar(model_params, [b, c, μ], gmm_pdf).zero[1]
end
μ_normal_vstars = map(μs) do μ
    FindVStar(model_params, [b, c, μ], normal_pdf).zero[1]
end

plot(
    μs, 
    μ_gmm_vstars, 
    label="Bimodal \$V^*(\\mu)\$", 
    xlabel="\$\\mu\$", 
    ylabel="\$V^*(\\mu)\$", 
    title="\$V^*\$ as \$\\mu\$ Changes")
plot!(
    μs, 
    μ_normal_vstars, 
    label="Gaussian \$V^*(\\mu)\$")

savefig(
    "temp-data/Vstar-mu.pdf"

)



plot(
    μs,
    1 .- gmm_cdf.(μ_gmm_vstars, Ref(model_params)),
    label = "Bimodal",
    xlabel = "\$\\mu\$",
    ylabel = "Take-up",
    title = "Take-up as \$\\mu\$ Changes"
)

plot!(
    μs,
    1 .- normal_cdf.(μ_normal_vstars, Ref(model_params)),
    label = "Normal",
    xlabel = "\$\\mu\$",
    ylabel = "Take-up"
)

savefig("temp-data/mu-changes.pdf")



function SocialMultiplier(params, fp_vector, pdf)
    b, c, μ = fp_vector
    vstar = FindVStar(params, fp_vector, pdf).zero[1]
    return 1 / (1 + μ * DeltaPrime(vstar, params, pdf))
end

plot(
    bs,
    map(bs) do b
        SocialMultiplier(model_params, [b, c, μ], gmm_pdf)
    end,
    label = "Bimodal",
    xlabel = "\$b\$",
    ylabel = "Social Multiplier",
    title = "Social Multiplier as \$b\$ Changes"
)


plot!(
    bs,
    map(bs) do b
        SocialMultiplier(model_params, [b, c, μ], normal_pdf)
    end,
    label = "Normal",
    xlabel = "\$b\$"
)
savefig("temp-data/social-multiplier-b.pdf")



small_μs = 0:0.01:0.5
plot(
    small_μs,
    map(small_μs) do μ
        SocialMultiplier(model_params, [b, c, μ], gmm_pdf)
    end,
    label = "Bimodal",
    xlabel = "\$\\mu\$",
    ylabel = "Social Multiplier",
    title = "Social Multiplier as \$\\mu\$ Changes"
)

plot!(
    small_μs,
    map(small_μs) do μ  
        SocialMultiplier(model_params, [b, c, μ], normal_pdf)
    end,
    label = "Normal",
    xlabel = "\$\\mu\$"
)
savefig("temp-data/social-multiplier-mu.pdf")



x = 0.5:0.01:2

plot(
    x,
    1 .- normal_cdf.(x, Ref(model_params)),
    label = "Normal Density for \$w\$"
)

x = -2:0.01:2
plot(
    x,
    1 .- gmm_cdf.(x, Ref(model_params)),
    label = "Bimodal Density for \$w\$"
)