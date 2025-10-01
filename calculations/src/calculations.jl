module calculations

using Unitful, PhysicalConstants.CODATA2018, CairoMakie, LaTeXStrings, Roots

function main()
    Tᴮ = 4501u"K"
    Tᴹ = 1789u"K"
    T_α_β = 244u"K"
    T_β_γ = 943u"K"
    ΔHⱽ = 121u"kJ/mol"
    ΔHᶠ = 65.3u"kJ/mol"
    ΔH_α_γ = 91.1u"kJ/mol"
    ΔH_β_α = 41.7u"kJ/mol"
    ΔH_α_β = -ΔH_β_α
    ΔH_β_γ = ΔH_α_γ - ΔH_β_α
    #  AW = 45.3u"g/mol"
    #  ρ_a = 5.45u"g/cm^3"
    #  ρ_β = 5.50u"g/cm^3"
    #  ρ_γ = 4.99u"g/cm^3"
    #  ρ_L = 4.28u"g/cm^3"
    n = 1
    R = BoltzmannConstant * AvogadroConstant

    #  ln_P_start = log(2e-58)
    ln_P_end = log(101325)

    # Liquid-Vapor
    #  LV_t_start = Tᴹ
    #  LV_t_end = Tᴮ
    LV_t_range = range(Tᴹ, Tᴮ, 100)
    LV_ln_P = @. -ΔHⱽ / (n * R) * (1 / LV_t_range - 1 / Tᴮ) + ln_P_end

    # Triple-point 1
    ln_P_min_1 = -ΔHⱽ / (n * R) * (1 / Tᴹ - 1 / Tᴮ) + ln_P_end

    # γ V coexistence
    ΔH_γ_V = ΔHᶠ + ΔHⱽ
    γL_t_range = range(T_β_γ, Tᴹ, 100)
    ln_P_γ_L = @. -ΔH_γ_V / (n * R) * (1 / γL_t_range - 1 / Tᴹ) + ln_P_min_1

    # β γ V triple point
    ln_P_min_2 = -ΔH_γ_V / (n * R) * (1 / T_β_γ - 1 / Tᴹ) + ln_P_min_1

    # γ V coexistence
    ΔH_β_V = ΔH_β_γ + ΔH_γ_V
    βγ_t_range = range(T_α_β, T_β_γ, 100)
    ln_P_βγ = @. -ΔH_β_V / (n * R) * (1 / βγ_t_range - 1 / T_β_γ) + ln_P_min_2

    # α β V triple point
    ln_P_min_3 = -ΔH_β_V / (n * R) * (1 / T_α_β - 1 / T_β_γ) + ln_P_min_2

    # α V coexistence
    ΔH_α_V = ΔH_α_γ + ΔHⱽ
    αβ_t_range = range(T_α_β / 1.2, T_α_β, 100)
    ln_P_αβ = @. -ΔH_α_V / (n * R) * (1 / αβ_t_range - 1 / T_α_β) + ln_P_min_3

    function diff_eqs(T)
        lhs = -(ΔH_γ_V / (n * R)) * (1 / (T * 1u"K") - 1 / Tᴹ)
        rhs = -(ΔH_α_V / (n * R)) * (1 / (T * 1u"K") - 1 / T_α_β)
        return lhs - rhs
    end


    T_solution = find_zero(diff_eqs, (1, 2000.0), Bisection())
    ln_P = -(ΔH_γ_V / (n * R)) * (1 / (T_solution * 1u"K") - 1 / Tᴹ) + log(ln_P_min_1)
    P_solution = exp(ln_P)
    println(T_solution)
    println(P_solution)

    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"1/x", ylabel=L"\ln(P)")
        ax.xreversed = true
        lv = lines!(ax, ustrip.(1 ./ LV_t_range), LV_ln_P, color=:red)
        tp_γlv = lines!(ax, [ustrip.(1 / Tᴹ), ustrip.(1 / Tᴹ)], [ln_P_min_1, ln_P_end], color=:blue)
        γl = lines!(ax, ustrip.(1 ./ γL_t_range), ln_P_γ_L, color=:green)
        tp_βγv = lines!(ax, [ustrip.(1 / T_β_γ), ustrip.(1 / T_β_γ)], [ln_P_min_2, ln_P_end], color=:orange)
        βγ = lines!(ax, ustrip.(1 ./ βγ_t_range), ln_P_βγ, color=:purple)
        tp_αβv = lines!(ax, [ustrip.(1 / T_α_β), ustrip.(1 / T_α_β)], [ln_P_min_3, ln_P_end], color=:magenta)
        αβ = lines!(ax, ustrip.(1 ./ αβ_t_range), ln_P_αβ, color=:maroon)

        text!(ax, 1 / Tᴹ * 2.6, -20,
            text="($(round(ustrip(Tᴹ), sigdigits=5)) K, $(round(exp(ln_P_min_1), sigdigits=5)) Pa)",
            align=(:left, :top), offset=(5, 5))
        text!(ax, (1 / Tᴹ) * 0.8, -40,
            text="V",
            align=(:right, :top), offset=(5, 5))
        text!(ax, (1 / Tᴹ) * 0.8, ln_P_min_1 * 1.5,
            text="L",
            align=(:right, :top), offset=(5, 5))
        text!(ax, (1 / Tᴹ) * 1.3, ln_P_min_1 * 1,
            text="γ",
            align=(:right, :top), offset=(5, 5))
        text!(ax, (1 / T_β_γ) * 1.1, ln_P_min_2 * 0.7,
            text="β",
            align=(:right, :top), offset=(5, 5))
        text!(ax, 1 / T_β_γ * 1.6, -30,
            text="($(round(ustrip(T_β_γ), sigdigits=5)) K, $(round(exp(ln_P_min_2), sigdigits=5)) Pa)",
            align=(:left, :top), offset=(5, 5))
        text!(ax, 1 / T_α_β, ln_P_min_3,
            text="($(round(ustrip(T_α_β), sigdigits=5)) K, $(round(exp(ln_P_min_3), sigdigits=5)) Pa)",
            align=(:left, :top), offset=(5, 5))
        text!(ax, (1 / T_α_β) * 1.1, ln_P_min_3 * 0.9,
            text="α",
            align=(:right, :top), offset=(5, 5))
        save("../assets/q2_fig.png", fig)
        fig
    end
end

end
