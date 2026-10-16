################################################################################
# Валидация  :luldt.
################################################################################

# --- независимая реализация правил -------------------------------------------
_lintrap(t₁, t₂, c₁, c₂) = (t₂ - t₁) * (c₁ + c₂) / 2
_logtrap(t₁, t₂, c₁, c₂) = (t₂ - t₁) * (c₂ - c₁) / log(c₂ / c₁)

"""
    trapauc(time, obs, rule) -> Float64

Прямая реализация правил трапеций по формулам, без обращения к внутренностям
MetidaNCA. Используется как независимый эталон.
"""
function trapauc(time, obs, rule)
    tmaxn = argmax(obs)                       # первый максимум, как в ctmax
    s = zero(eltype(obs)) * zero(eltype(time))
    for i in 1:length(time) - 1
        t₁, t₂, c₁, c₂ = time[i], time[i+1], obs[i], obs[i+1]
        after = i >= tmaxn
        uselog =
            rule === :lint  ? false :
            rule === :logt  ? (after && c₁ > 0 && c₂ > 0) :
            rule === :luld  ? (c₁ > c₂ > 0) :
            rule === :luldt ? (after && c₁ > c₂ > 0) :
            error("Unknown rule: $rule")
        s += uselog ? _logtrap(t₁, t₂, c₁, c₂) : _lintrap(t₁, t₂, c₁, c₂)
    end
    return s
end

@testset verbose = true "  Trapezoidal rules; analytic validation                  " begin

    # ---------------------------------------------------------------- T1 ----
    # Лог-трапеция на экспоненте точна:
    #   (t₂-t₁)(C₂-C₁)/ln(C₂/C₁) = (C₁-C₂)/k = ∫C dt
    @testset "  T1 Log rules are exact on a mono-exponential profile    " begin
        k, C0 = 0.35, 100.0
        t = [0.0, 0.5, 1.0, 2.0, 3.5, 5.0, 8.0, 12.0]
        c = C0 .* exp.(-k .* t)
        analytic = (first(c) - last(c)) / k

        for rule in (:logt, :luld, :luldt)
            @test MetidaNCA.auc(t, c; calcm = rule) ≈ analytic rtol = 1e-12
        end
        # линейная трапеция на выпуклой убывающей функции завышает площадь
        @test MetidaNCA.auc(t, c; calcm = :lint) > analytic
    end

    # ---------------------------------------------------------------- T2 ----
    # Та же экспонента через полный расчёт nca!
    # AUClast, AUMClast, AUCinf, AUMCinf, Kel, HL, MRTinf.
    @testset "  T2 Full NCA chain against closed-form values            " begin
        k, C0 = 0.35, 100.0
        t = [0.0, 0.5, 1.0, 2.0, 3.5, 5.0, 8.0, 12.0]
        c = C0 .* exp.(-k .* t)
        tl, cl = last(t), last(c)

        exp_auclast  = (C0 - cl) / k
        exp_aumclast = -tl * cl / k + (C0 - cl) / k^2      # ∫t·C dt от 0 до tl
        exp_aucinf   = C0 / k                              # ∫ от 0 до ∞
        exp_aumcinf  = C0 / k^2

        for rule in (:logt, :luld, :luldt)
            subj = MetidaNCA.pkimport(t, c;
                       dosetime = MetidaNCA.DoseTime(dose = 1.0, time = 0.0))
            r = MetidaNCA.nca!(subj; adm = :iv, calcm = rule)

            @test r[:Rsq]       ≈ 1.0            rtol = 1e-12   # точная прямая в логарифмах
            @test r[:Kel]       ≈ k              rtol = 1e-12
            @test r[:HL]        ≈ log(2) / k     rtol = 1e-12
            @test r[:AUClast]   ≈ exp_auclast    rtol = 1e-12
            @test r[:AUMClast]  ≈ exp_aumclast   rtol = 1e-12
            @test r[:AUCinf]    ≈ exp_aucinf     rtol = 1e-12
            @test r[:AUMCinf]   ≈ exp_aumcinf    rtol = 1e-12
            @test r[:MRTinf]    ≈ 1 / k          rtol = 1e-12
        end
    end

    # ---------------------------------------------------------------- T3 ----
    # Строго возрастающий профиль: лог не применяется ни одним правилом,
    # все четыре обязаны дать одно и то же - площадь линейных трапеций.
    @testset "  T3 Linear branch on a strictly ascending profile        " begin
        t = [0.0, 1.0, 2.5, 4.0, 7.0]
        c = [0.0, 20.0, 55.0, 80.0, 140.0]
        lin = sum(_lintrap(t[i], t[i+1], c[i], c[i+1]) for i in 1:length(t)-1)

        for rule in (:lint, :logt, :luld, :luldt)
            @test MetidaNCA.auc(t, c; calcm = rule) ≈ lin rtol = 1e-12
        end
    end

    # ---------------------------------------------------------------- T4 ----
    # Правило переключения через тождества с внешне валидированными методами.
    @testset "  T4 Switching-rule identities                            " begin
        auc(t, c, m) = MetidaNCA.auc(t, c; calcm = m)

        # A. один пик: рост до Tmax, затем спад.
        #    До Tmax лог не нужен, после Tmax спад - все три берут лог.
        tA = [0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
        cA = [0.0, 50.0, 100.0, 70.0, 40.0, 15.0, 5.0]
        @test auc(tA, cA, :luldt) ≈ auc(tA, cA, :logt) rtol = 1e-12
        @test auc(tA, cA, :luldt) ≈ auc(tA, cA, :luld) rtol = 1e-12

        # B. двугорбый: есть ВОЗРАСТАЮЩИЙ участок после Tmax.
        #    Единственный случай, где :luldt расходится с :logt.
        tB = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0]
        cB = [0.0, 50.0, 100.0, 70.0, 40.0, 60.0, 45.0, 20.0, 6.0]
        @test !isapprox(auc(tB, cB, :luldt), auc(tB, cB, :logt); rtol = 1e-12)
        @test auc(tB, cB, :luldt) ≈ auc(tB, cB, :luld) rtol = 1e-12

        # C. монотонный спад с первой точки: Tmax - первая точка, всё после неё.
        tC = [0.0, 1.0, 2.0, 4.0, 6.0, 10.0]
        cC = [100.0, 70.0, 50.0, 25.0, 12.0, 4.0]
        @test auc(tC, cC, :luldt) ≈ auc(tC, cC, :logt) rtol = 1e-12
        @test auc(tC, cC, :luldt) ≈ auc(tC, cC, :luld) rtol = 1e-12

        # D. провал ДО Tmax.
        #    Единственный случай, где :luldt расходится с :luld.
        tD = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 9.0, 14.0]
        cD = [0.0, 80.0, 45.0, 120.0, 70.0, 40.0, 18.0, 6.0]
        @test !isapprox(auc(tD, cD, :luldt), auc(tD, cD, :luld); rtol = 1e-12)
        @test auc(tD, cD, :luldt) ≈ auc(tD, cD, :logt) rtol = 1e-12

        # контроль: профили действительно содержат нужные участки
        @test any(i -> i >= argmax(cB) && cB[i] < cB[i+1], 1:length(cB)-1)   # рост после Tmax
        @test any(i -> i <  argmax(cD) && cD[i] > cD[i+1], 1:length(cD)-1)   # спад до Tmax
    end

    # ---------------------------------------------------------------- T5 ----
    # Независимый пересчёт на реальных данных
    @testset "  T5 Independent reimplementation on pkdata2              " begin
        for subj in sort(unique(pkdata2.Subject))
            g = pkdata2[pkdata2.Subject .== subj, :]
            t = float.(g.Time)
            c = float.(g.Concentration)
            for rule in (:lint, :logt, :luld, :luldt)
                @test MetidaNCA.auc(t, c; calcm = rule) ≈ trapauc(t, c, rule) rtol = 1e-12
            end
        end
    end
end
