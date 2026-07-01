using DifferentialEquations
using Plots
using Statistics

# ── Model ─────────────────────────────────────────────────────────────────────
function dde_model!(du, u, h, p, t)
    alpha, tau, gamma, c = p.alpha, p.tau, p.gamma, p.c
    x         = u[1]
    x_delayed = h(p, t - tau)[1]
    du[1] = x - x^3 - alpha * x_delayed * (1 - gamma * x_delayed^2) + c * cos(2π * t)
end

# ── Find local maxima directly on solver output points ────────────────────────
# Much faster than building a 50k uniform grid.
function find_local_maxima(x::Vector{Float64})
    maxima = Float64[]
    @inbounds for i in 2:length(x)-1
        if x[i] > x[i-1] && x[i] > x[i+1]
            push!(maxima, x[i])
        end
    end
    return maxima
end

# ── Bifurcation diagram ───────────────────────────────────────────────────────
function bifurcation_diagram()
    tau_values     = range(4.1, 3.91, length=3000)
    transient_time = 1500.0          # was 6000 — DDE attractors settle faster
    tspan          = (0.0, 5000.0)   # was 10000 — only integrate what you use

    tau_plot = Float64[]
    max_plot = Float64[]

    # Constant history — no prev_sol trick (see notes)
    history_func = (p, t) -> [1.3]

    for tau in tau_values
        print("τ = $(round(tau, digits=4)) \r")

        params = (alpha=0.75, tau=tau, gamma=0.49, c=1.52, x0=1.3)

        prob = DDEProblem(
            dde_model!,
            [params.x0],
            history_func,
            tspan,
            params;
            constant_lags = [tau]
        )

        # reltol 1e-6 is plenty for a bifurcation diagram
        sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-8, abstol=1e-8)

        # ── Use the solver's own adaptive output points ────────────────────
        # Filter to post-transient, then find maxima — no uniform grid needed.
        idx     = findall(>=(transient_time), sol.t)
        x_post  = [sol.u[i][1] for i in idx]
        maxima  = find_local_maxima(x_post)

        if !isempty(maxima)
            append!(tau_plot, fill(tau, length(maxima)))
            append!(max_plot, maxima)
        end
    end

    println("\nDone!")

    fig = scatter(tau_plot, max_plot;
                  markersize   = 0.0005,
                  markerstroke = stroke(0),
                  markeralpha   = 0.045,  # 0 = fully transparent, 1 = opaque
                  label        = "",
                  xlabel       = "τ",
                  ylabel       = "local maxima of x(t)",
                  title        = "Bifurcation diagram  (α=0.75, γ=0.49, c=1.52)",
                  dpi          = 300,
                  grid         = false)
                  #xticks       = false,      # remove x-axis ticks
                  #yticks       = false,      # remove y-axis ticks
                  #xaxis        = false,      # remove x-axis line and ticks (optional)
                  #yaxis        = false)
    display(fig)
    savefig(fig, "chaosinMZ5.png")

    return tau_plot, max_plot
end

tau_plot, max_plot = bifurcation_diagram()